function results = runConvergence51(varargin)
% convergence test 5.1 of Keilegavlen & Nordbotten (2017)
%
%   runConvergence51('dim',3,'grid','simplex','nu',[0.25 0.495])
%   runConvergence51('dim', 2, 'grid', 'cartesian', 'material', 'heterogeneous', 'kappa', 1e6, 'nu', 0.25, 'levels', [2 4 8 16 32]);
%
% options
%   dim      2 or 3                                  (default 2)
%   grid     'cartesian' | 'simplex'                 (default 'cartesian')
%   levels   vector of cells per direction           (default [2 4 8 16])
%   material 'heterogeneous' (default) keeps the corner contrast of section
%            5.1.1 with the given kappa; 'homogeneous' removes it entirely by
%            forcing kappa = 1, so the ONLY thing that varies over the sweep is
%            the Poisson ratio
%   kappa    heterogeneity contrast, used only when material is
%            'heterogeneous'                                (default 1e6)
%   nu       Poisson ratio, SCALAR OR VECTOR. one convergence study is run per
%            value, with  alpha = lambda/mu = 2 nu / (1 - 2 nu)
%            (default [0.25 0.495 0.4995 0.49995], i.e. alpha ~ 1, 1e2, 1e3, 1e4,
%            the incompressible sweep of section 5.1.4 plus a reference)
%   alpha    used only when nu is empty, for a direct lambda/mu
%   perturb  interior vertex perturbation amplitude  (default 0)
%   seed     rng seed for the perturbation           (default 0)
%   solver   'iterative' (default), 'direct' or 'gmres'.
%              iterative : MINRES with a block incomplete-Cholesky
%                          preconditioner, O(n) memory (solveSaddleLean)
%              direct    : backslash; pays the fill-in of the saddle system,
%                          which is expensive from 3D N = 16 (4096 cells) on
%              gmres     : solvers/iterative_solver, right-preconditioned GMRES
%                          with exact factorizations; the heaviest of the three
%   quadDegree  exactness degree of the cell/face quadrature (default 6).
%            the 5.1 solution is degree 6, so the Dirichlet integrand is
%            degree 7; the built-in degree-2 rules under-integrate it
%   problem  supply a problem struct to override 5.1 entirely (nu is ignored)
%   bc       'auto' (default), 'dirichlet' or 'slab'. 'auto' uses 'slab' for
%            dim = 2 and 'dirichlet' for dim = 3. 'slab' gives the artificial
%            z-normal faces of the quasi-2D slab the exact traction as an
%            essential condition instead of leaving sigma_zz free, which would
%            otherwise absorb lambda*tr(eps_h) = O(lambda h^2)
%   verbose  print the tables                        (default true)
%
% dim = 2 solves the plane problem on a one-cell-thick slab: the 2D solution
% extended by u_z = 0 with no z dependence is an exact 3D solution, so the
% existing 3D discretization is used unchanged.
%
% errors are the four norms of the paper, computed by computeErrorsDassi:
%
%   E_u      = || u - u_h ||_0
%   E_sigdiv = || div sigma - div sigma_h ||_0
%   E_sigPi  = || sigma - Pi_E sigma_h ||_0
%   E_sig    = ( sum_f h_f int_f kappa |(sigma - sigma_h) n|^2 df )^1/2
%
% results is a struct array, one entry per nu, each holding the per-level
% vectors N, h, nCells, nDofs, eU, eDiv, eProj, eSig, time.

    setupVerificationPath();

    opt = struct('dim',2,'grid','cartesian','levels',[2 4 8 16], ...
                 'material','heterogeneous','kappa',1e6, ...
                 'alpha',1,'nu',[0.25 0.495 0.4995 0.49995], ...
                 'perturb',0,'seed',0, ...
                 'solver','iterative','verbose',true,'problem',[],'quadDegree',6, ...
                 'bc','auto');

    for k = 1:2:numel(varargin)
        assert(isfield(opt, varargin{k}), 'Unknown option "%s".', varargin{k});
        opt.(varargin{k}) = varargin{k+1};
    end

    if strcmp(opt.bc, 'auto')
        if opt.dim == 2 && isempty(opt.problem)
            opt.bc = 'slab';
        else
            opt.bc = 'dirichlet';
        end
    end

    switch opt.material
        case 'homogeneous'
            opt.kappa = 1;             % chi cancels: mu = 1, lambda = alpha
        case 'heterogeneous'
            % keep opt.kappa as given
        otherwise
            error('Unknown material "%s". Use ''homogeneous'' or ''heterogeneous''.', opt.material);
    end

    % a supplied problem fixes the material, so nu plays no role
    if ~isempty(opt.problem)
        nuList = NaN;
    elseif ~isempty(opt.nu)
        nuList = opt.nu(:).';
    else
        nuList = NaN;                  % fall back to the explicit alpha
    end

    results = repmat(struct('nu',[],'alpha',[],'N',[],'h',[],'nCells',[], ...
                            'nDofs',[],'eU',[],'eDiv',[],'eProj',[],'eSig',[], ...
                            'time',[]), 1, numel(nuList));

    for iNu = 1:numel(nuList)

        nu = nuList(iNu);

        if isnan(nu)
            alpha = opt.alpha;
        else
            alpha = 2*nu / (1 - 2*nu);
        end

        if isempty(opt.problem)
            problem = problemKeilegavlen51(opt.dim, opt.kappa, alpha);
        else
            problem = opt.problem;
        end

        r = runOneStudy(opt, problem, nu, alpha);

        results(iNu) = r;

    end

    if opt.verbose && numel(nuList) > 1
        printSummary(results);
    end

    plotConvergence51(results);

end

%% ------------------------------------------------------------------

function r = runOneStudy(opt, problem, nu, alpha)
% one full refinement study at a fixed material

    nL = numel(opt.levels);

    r.nu = nu;
    r.alpha = alpha;
    r.N = opt.levels(:).';
    r.h = 1 ./ r.N;
    r.nCells = zeros(1,nL);
    r.nDofs = zeros(1,nL);
    r.eU = zeros(1,nL);
    r.eDiv = zeros(1,nL);
    r.eProj = zeros(1,nL);
    r.eSig = zeros(1,nL);
    r.time = zeros(1,nL);

    if opt.verbose
        if isnan(nu)
            head = problem.name;
        else
            head = sprintf('%s | nu = %g', problem.name, nu);
        end
        fprintf('\n%s | %s grid | perturb %.3g | bc %s | solver %s\n', ...
                head, opt.grid, opt.perturb, opt.bc, opt.solver);
        fprintf('%9s %7s %7s %11s %6s %11s %6s %11s %6s %11s %6s\n', ...
                'dofs','cells','h','E_u','rate','E_sig,div','rate', ...
                'E_sig,Pi','rate','E_sig','rate');
    end

    for k = 1:nL

        t0 = tic;

        [cell_struct, face_struct, V3, cells3D] = makeMesh(opt, r.N(k));

        cell_struct = assignCellPropertiesFun(cell_struct, V3, cells3D, problem);

        switch opt.bc
            case 'dirichlet'
                face_struct = assignFacePropertiesDirichlet(face_struct);
            case 'slab'
                face_struct = assignFacePropertiesSlab(face_struct);
            otherwise
                error('Unknown bc "%s".', opt.bc);
        end

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

        % only B_local, D_local and P_local are needed after the solve, for the
        % error norms. releasing the rest before the factorization keeps the
        % peak well below the direct solver's own fill-in
        clear A_local Bproj K_cons K_stab

        nS = 6*numel(face_struct);
        nD = 6*numel(cell_struct);

        A = [K_global, B_global.'; B_global, sparse(nD,nD)];

        clear K_global B_global

        rhs = [assembleDirichletRHSFun(cell_struct, face_struct, face_global_geom, problem);
               assembleBodyForce(cell_struct, problem)];

        if any([face_struct.is_neumann])
            [A, rhs] = assembleNeumannFun(face_struct, A, rhs, face_global_geom, problem);
        end

        switch opt.solver
            case 'direct'
                sol = A \ rhs;
            case 'iterative'
                [sol, sinfo] = solveSaddleLean(A, rhs, nS, nD);
                if sinfo.relres > 1e-8
                    warning('N=%d: iterative solver reached relres %.2e only (flag %d, %d it)', ...
                            r.N(k), sinfo.relres, sinfo.flag, sinfo.iterations);
                end
            case 'gmres'
                % solvers/iterative_solver: right-preconditioned GMRES with exact
                % factorizations of K and of the approximate Schur complement.
                % much heavier than 'iterative' in both memory and time
                [s_h, u_hh] = iterative_solver(A, rhs, cell_struct, face_struct);
                sol = [s_h; u_hh];
            otherwise
                error('Unknown solver "%s".', opt.solver);
        end

        clear A

        err = computeErrorsDassi(cell_struct, face_struct, face_global_geom, V3, ...
                                 B_local, D_local, P_local, ...
                                 sol(1:nS), sol(nS+1:end), problem);

        r.nCells(k) = numel(cell_struct);
        r.nDofs(k) = nS + nD;
        r.eU(k) = err.displacement;
        r.eDiv(k) = err.divergence;
        r.eProj(k) = err.projection;
        r.eSig(k) = err.stress;
        r.time(k) = toc(t0);

        if opt.verbose
            fprintf('%9d %7d %7.4f %11.4e %6s %11.4e %6s %11.4e %6s %11.4e %6s\n', ...
                r.nDofs(k), r.nCells(k), r.h(k), ...
                r.eU(k),    rateStr(r.h, r.eU,    k), ...
                r.eDiv(k),  rateStr(r.h, r.eDiv,  k), ...
                r.eProj(k), rateStr(r.h, r.eProj, k), ...
                r.eSig(k),  rateStr(r.h, r.eSig,  k));
        end

    end

end

function s = rateStr(h, e, k)
% observed rate between levels k-1 and k

    if k == 1 || ~isfinite(e(k)) || ~isfinite(e(k-1)) || e(k) <= 0 || e(k-1) <= 0
        s = '  -  ';
    else
        s = sprintf('%5.2f', log(e(k-1)/e(k)) / log(h(k-1)/h(k)));
    end

end

function printSummary(results)
% asymptotic rates over the finest two levels, one row per nu

    fprintf('\n=== summary: rate over the finest two levels ===\n');
    fprintf('%9s %10s %9s %8s %8s %8s %8s %13s %13s\n', ...
            'nu','alpha','dofs','r_u','r_sdiv','r_sPi','r_sig','E_u','E_sig');

    for i = 1:numel(results)

        r = results(i);
        k = numel(r.N);

        rate = @(e) log(e(k-1)/e(k)) / log(r.h(k-1)/r.h(k));

        if k < 2
            ru = NaN; rd = NaN; rp = NaN; rs = NaN;
        else
            ru = rate(r.eU); rd = rate(r.eDiv); rp = rate(r.eProj); rs = rate(r.eSig);
        end

        fprintf('%9.5f %10.4g %9d %8.2f %8.2f %8.2f %8.2f %13.4e %13.4e\n', ...
                r.nu, r.alpha, r.nDofs(k), ru, rd, rp, rs, r.eU(k), r.eSig(k));

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

    needed = {'geometry','operators','assembly','solvers','io','meshes','verification_incompressibility'};

    for k = 1:numel(needed)
        d = fullfile(root, needed{k});
        if isfolder(d) && ~contains([path pathsep], [d pathsep])
            addpath(d);
        end
    end

end

function plotConvergence51(results)

    for i = 1:numel(results)

        r = results(i);

        figure('Color', 'w');

        loglog(r.h, r.eU,'-o', 'LineWidth', 1.5); 
        hold on;
        loglog(r.h, r.eDiv, '-s', 'LineWidth', 1.5);
        loglog(r.h, r.eProj, '-^', 'LineWidth', 1.5);
        loglog(r.h, r.eSig, '-d', 'LineWidth', 1.5);

        % O(h) reference
        ref = r.eSig(end) / r.h(end) * r.h;
        loglog(r.h, ref, '--', 'LineWidth', 1.2);

        grid on;
        box on;

        xlim([min(r.h)/1.1, max(r.h)*1.1]);
        xticks(sort(r.h));
        xticklabels(compose('%.4g', sort(r.h)));

        xlabel('$h$', 'Interpreter', 'latex');
        ylabel('Error');

        legend('$E_u$', '$E_{\sigma,\mathrm{div}}$', '$E_{\sigma,\Pi}$', ...
               '$E_\sigma$', ...
               '$O(h)$', ...
               'Interpreter', 'latex', ...
               'Location', 'best');

        if ~isempty(r.nu) && ~isnan(r.nu)
            title(sprintf('Convergence Study ($\\nu = %.5g$)', r.nu), 'Interpreter', 'latex');
        else
            title('Convergence Study');
        end

        set(gca, 'FontSize', 12, 'TickLabelInterpreter', 'latex');

        hold off;

    end

end