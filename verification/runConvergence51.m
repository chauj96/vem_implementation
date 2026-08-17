function results = runConvergence51(varargin)
% convergence test 5.1 of Keilegavlen & Nordbotten (2017)
%
%   results = runConvergence51('dim',3,'grid','cartesian','levels',[2 4 8])
%
% options
%   dim      2 or 3                                  (default 3)
%   grid     'cartesian' | 'simplex'                 (default 'cartesian')
%   levels   vector of cells per direction           (default [2 4 8 16])
%   kappa    heterogeneity contrast                  (default 1)
%   alpha    lambda / mu                             (default 1)
%   nu       Poisson ratio; if given, overrides alpha via
%            alpha = 2 nu / (1 - 2 nu)   (nu = 0.25 <-> alpha = 1)
%   perturb  interior vertex perturbation amplitude  (default 0)
%   seed     rng seed for the perturbation           (default 0)
%   solver   'direct' | 'iterative'                  (default 'direct')
%   quadDegree  exactness degree of the cell/face quadrature (default 6).
%            the 5.1 solution is degree 6, so the Dirichlet integrand is
%            degree 7; the built-in degree-2 rules under-integrate it
%   verbose  print the table                         (default true)
%
% dim = 2 solves the plane problem on a one-cell-thick slab: the 2D solution
% extended by u_z = 0 with no z dependence is an exact 3D solution, so the
% existing 3D discretization is used unchanged.
%
% this driver does not modify any existing routine. it swaps in the
% heterogeneous material, the analytical Dirichlet data and the body force,
% all of which live in verification/

    setupVerificationPath();

    opt = struct('dim',2,'grid','cartesian','levels',[2 4 8 16], ...
                 'kappa',1,'alpha',1,'nu',[0.495],'perturb',0,'seed',0, ...
                 'solver','direct','verbose',true,'problem',[],'quadDegree',6);

    for k = 1:2:numel(varargin)
        assert(isfield(opt, varargin{k}), 'Unknown option "%s".', varargin{k});
        opt.(varargin{k}) = varargin{k+1};
    end

    if ~isempty(opt.nu)
        opt.alpha = 2*opt.nu / (1 - 2*opt.nu);
    end

    if isempty(opt.problem)
        problem = problemKeilegavlen51(opt.dim, opt.kappa, opt.alpha);
    else
        problem = opt.problem;      % e.g. problemPatchTest, for validating the driver
    end

    nL = numel(opt.levels);
    results = struct('N',cell(nL,1),'h',[],'nCells',[],'nDofs',[], ...
                     'eD',[],'eS',[],'eR',[],'time',[]);

    if opt.verbose
        fprintf('\n%s | %s grid | perturb %.3g\n', problem.name, opt.grid, opt.perturb);
        fprintf('%5s %8s %8s %12s %7s %12s %7s %12s %7s\n', ...
                'N','cells','h','e_disp','rate','e_stress','rate','e_rot','rate');
    end

    for k = 1:nL

        N = opt.levels(k);
        t0 = tic;

        [cell_struct, face_struct, V3, cells3D] = makeMesh(opt, N);

        cell_struct = assignCellPropertiesFun(cell_struct, V3, cells3D, problem);
        face_struct = assignFacePropertiesDirichlet(face_struct);

        cell_struct = computeVolumeQuadratureP(cell_struct, face_struct, V3, opt.quadDegree);
        face_struct = computeFaceQuadratureP(face_struct, V3, opt.quadDegree);

        A_local = computeAE_Tonon(cell_struct, face_struct, V3);
        B_local = computeBE(cell_struct, face_struct, V3);
        face_global_geom = computeGlobalFaceGeometry(face_struct, V3);
        D_local = computeDE(cell_struct, A_local, B_local);
        Bproj = computeBproj(cell_struct, face_struct, D_local, B_local);
        P_local = computeProjectionMatrix(cell_struct, Bproj);
        K_cons = computeConsistency(cell_struct, P_local);
        K_stab = computeStabilization(cell_struct, face_struct, B_local, P_local);

        [K_global, B_global] = assembleSaddleSystem(cell_struct, face_struct, ...
                                   face_global_geom, B_local, K_cons, K_stab);

        nS = 6*numel(face_struct);
        nD = 6*numel(cell_struct);

        A = [K_global, B_global.'; B_global, sparse(nD,nD)];

        rhs = [assembleDirichletRHSFun(cell_struct, face_struct, face_global_geom, problem);
               assembleBodyForce(cell_struct, problem)];

        switch opt.solver
            case 'direct'
                sol = A \ rhs;
            case 'iterative'
                [s_h, u_hh] = iterative_solver(A, rhs, cell_struct, face_struct);
                sol = [s_h; u_hh];
            otherwise
                error('Unknown solver "%s".', opt.solver);
        end

        sigma_h = sol(1:nS);
        u_h = sol(nS+1:end);

        err = computeErrors51(cell_struct, face_struct, face_global_geom, sigma_h, u_h, problem);

        results(k).N = N;
        results(k).h = 1/N;
        results(k).nCells = numel(cell_struct);
        results(k).nDofs = nS + nD;
        results(k).eD = err.displacement;
        results(k).eS = err.stress;
        results(k).eR = err.rotation;
        results(k).time = toc(t0);

        if opt.verbose
            if k == 1
                rD = NaN; rS = NaN; rR = NaN;
            else
                lr = log(results(k-1).h / results(k).h);
                rD = log(results(k-1).eD / results(k).eD) / lr;
                rS = log(results(k-1).eS / results(k).eS) / lr;
                rR = log(results(k-1).eR / results(k).eR) / lr;
            end
            fprintf('%5d %8d %8.4f %12.4e %7.2f %12.4e %7.2f %12.4e %7.2f\n', ...
                N, results(k).nCells, results(k).h, ...
                results(k).eD, rD, results(k).eS, rS, results(k).eR, rR);
        end

    end

end

function [cell_struct, face_struct, V3, cells3D] = makeMesh(opt, N)
% dim = 2 uses a single cell layer in z, with thickness h so cells stay isotropic

    if opt.dim == 2
        n = [N N 1];
        L = [1 1 1/N];
    else
        n = [N N N];
        L = [1 1 1];
    end

    switch opt.grid
        case 'cartesian'
            [cell_struct, face_struct, V3, cells3D] = cartesianMeshBox(n, L, opt.perturb, opt.seed);
        case 'simplex'
            [cell_struct, face_struct, V3, cells3D] = simplexMeshBox(n, L, opt.perturb, opt.seed);
        otherwise
            error('Unknown grid "%s".', opt.grid);
    end

end

function setupVerificationPath()
% add the repository folders this driver needs, so it runs from a bare session

    root = fileparts(fileparts(mfilename('fullpath')));

    needed = {'geometry','operators','assembly','solvers','io','meshes','verification'};

    for k = 1:numel(needed)
        d = fullfile(root, needed{k});
        if isfolder(d) && ~contains([path pathsep], [d pathsep])
            addpath(d);
        end
    end

end
