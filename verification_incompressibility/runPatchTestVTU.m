function [err, results] = runPatchTestVTU(vtuFile, varargin)
% linear-displacement patch test on a mesh read from a .vtu file
%
%   runPatchTestVTU('meshes/polyhedral_voronoi/polyhedral_voronoi_complex.vtu')
%
% A single boundary-value problem. The exact displacement is the linear field
%
%   u(x) = G (x - x0) / L
%
% with x0 the centre of the mesh bounding box and L its characteristic length,
% so u is O(1) over the domain whatever the mesh units.
%
% G defaults to the identity: a uniform isotropic dilation, u = (x - x0)/L.
% eps = I/L and sigma = (2 mu + 3 lambda) I / L are constant, and the rotation
% is zero, so u_h reduces to a pure translation on each element. pass any other
% 3x3 G to exercise shear and rotation as well.
%
% grad u = G/L is constant, hence eps = sym(G)/L and sigma = 2 mu eps +
% lambda tr(eps) I are constant and f = div sigma = 0.
%
% Boundary conditions are pure Dirichlet on the whole external boundary. The
% external facets are found from the face -> cell incidence graph: a facet
% carried by a single cell is external. u is prescribed there.
%
% sigma is constant, so it lies in the stress space exactly. u is linear, which
% is NOT in RM(E), so the reachable target for u_h is the L2 projection Pi_RM u
% and E_u is measured against that. Every norm reported is therefore a
% consistency check and sits at the floor set by the linear solve.
%
% Dassi's ||u - u_h||_0 is reported too, but only for information: it measures
% the distance from RM(E) to a linear field and is O(h) by construction.
%
% options
%   G          displacement gradient, 3x3   (default eye(3))
%   lambda, mu                              (default 1, 1)
%   solver     'iterative' | 'direct'       (default 'iterative')
%   solverTol  tolerance of the iterative solve  (default 1e-5)
%   quadDegree                              (default 2; the fields are linear)
%   tol        pass threshold on the error norms  (default 1e-4). the solve sets
%              the attainable floor: every norm scales linearly with solverTol,
%              and E_sig,div lands about 20x above it because it differentiates
%              the stress and so amplifies the residual
%   output     .vtu file to write, '' to skip  (default 'output/patch_test.vtu')
%   verbose                                 (default true)
%
% The output carries six cell-data fields: displacement and displacement_exact
% (3 components), stress and stress_exact (9 components, row-major), and
% displacement_error and stress_error, which are the local Dassi L2 norms
% ||Pi_RM u - u_h||_{0,E} and ||sigma - Pi_E sigma_h||_{0,E}.

    opt = struct('G', eye(3), 'lambda',1, 'mu',1, ...
                 'solver','iterative', 'solverTol',1e-5, 'quadDegree',2, ...
                 'tol',1e-4, 'output',fullfile('output','patch_test.vtu'), ...
                 'verbose',true);

    for k = 1:2:numel(varargin)
        assert(isfield(opt, varargin{k}), 'Unknown option "%s".', varargin{k});
        opt.(varargin{k}) = varargin{k+1};
    end

    setupPatchPath();

    %% mesh and topology
    ev = evalc('[cell_struct, face_struct, V3, cells3D] = readVTU(vtuFile);'); %#ok<NASGU>

    topo = faceCellTopology(cell_struct, face_struct);

    lo = min(V3, [], 1);
    hi = max(V3, [], 1);

    box.centre = (lo + hi) / 2;
    box.extent = hi - lo;
    box.length = max(box.extent);          % characteristic length
    box.diagonal = norm(box.extent);

    if opt.verbose
        fprintf('\n=== linear patch test: %s ===\n', vtuFile);
        fprintf('cells %d, faces %d, vertices %d\n', ...
                numel(cell_struct), numel(face_struct), size(V3,1));
        fprintf('face->cell graph: %d external (1 cell), %d internal (2 cells)\n', ...
                topo.nBoundary, topo.nInterior);
        fprintf('bounding box [%.4g %.4g %.4g] x [%.4g %.4g %.4g], L = %.4g\n', ...
                lo, hi, box.length);
    end

    assert(isempty(topo.orphan), ...
        '%d facets have no incident cell.', numel(topo.orphan));
    assert(isempty(topo.nonManifold), ...
        '%d facets have more than two incident cells; mesh is non-manifold.', numel(topo.nonManifold));
    assert(isempty(topo.mismatch), ...
        '%d facets disagree between the cell->face lists and face_struct.cells.', numel(topo.mismatch));
    assert(topo.nBoundary > 0, 'No external facets found.');

    %% solve
    problem = problemLinearPatch(opt.G, box, opt.lambda, opt.mu, 'linear patch');

    [err, cellData, cellOut] = solveOne(cell_struct, face_struct, V3, cells3D, problem, opt);

    if ~isempty(opt.output)
        writeMeshVTU(opt.output, V3, cellOut, face_struct, cellData);
    end

    results = struct('G', opt.G, 'box', box, 'topology', topo, ...
                     'err', err, 'cellData', cellData);
    results.pass = max([err.displacementRM, err.divergence, ...
                        err.projection, err.stress]) < opt.tol;

    if opt.verbose

        fprintf('\n%-14s %11s %11s %11s %11s\n', ...
                '', 'E_u(Pi_RM)', 'E_sig,div', 'E_sig,Pi', 'E_sig');
        fprintf('%-14s %11.3e %11.3e %11.3e %11.3e\n', 'error', ...
                err.displacementRM, err.divergence, err.projection, err.stress);

        fprintf('\n%s (all norms below %.0e, solverTol %.0e)\n', ...
                ternary(results.pass,'PASS','FAIL'), opt.tol, opt.solverTol);
        fprintf('||u - u_h||_0 = %.3e, Dassi''s norm against the exact field:\n', ...
                err.displacement);
        fprintf('O(h) by construction, not part of the verdict.\n');

    end

end

%% ------------------------------------------------------------------

function [err, cd, cell_struct] = solveOne(cell_struct, face_struct, V3, cells3D, problem, opt)

    cell_struct = assignCellPropertiesFun(cell_struct, V3, cells3D, problem);

    % pure Dirichlet: every facet with a single incident cell
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

    clear A_local Bproj K_cons K_stab

    nS = 6*numel(face_struct);
    nD = 6*numel(cell_struct);

    A = [K_global, B_global.'; B_global, sparse(nD,nD)];
    clear K_global B_global

    rhs = [assembleDirichletRHSFun(cell_struct, face_struct, face_global_geom, problem);
           assembleBodyForce(cell_struct, problem)];

    switch opt.solver
        case 'direct'
            sol = A \ rhs;
        case 'iterative'
            sol = solveSaddleLean(A, rhs, nS, nD, opt.solverTol);
        otherwise
            error('Unknown solver "%s".', opt.solver);
    end

    clear A

    sigma_h = sol(1:nS);
    u_h = sol(nS+1:end);

    err = computeErrorsDassi(cell_struct, face_struct, face_global_geom, V3, ...
                             B_local, D_local, P_local, sigma_h, u_h, problem);

    err.displacementRM = rmProjectionError(cell_struct, u_h, problem);

    cd = patchCellData(cell_struct, face_struct, face_global_geom, ...
                       B_local, P_local, sigma_h, u_h, problem);

end

function eRel = rmProjectionError(cell_struct, u_h, problem)
% || Pi_RM u - u_h ||_0 / || u ||_0
%
% RM(E) is spanned by phi_i = e_i and phi_{3+i} = e_i x (x - x_E). u_h has
% coefficients [alpha; omega] in exactly this basis, so the error reduces to
% (c - d)' M (c - d) with M the local RM Gram matrix and c the coefficients of
% the L2 projection of the exact field, M c = b.

    e2 = 0; n2 = 0;

    for e = 1:numel(cell_struct)

        qp = cell_struct(e).quad_points;
        qw = cell_struct(e).quad_weights(:).';
        r = qp - cell_struct(e).center(:);

        ue = problem.u(qp);

        [M, b] = rmGram(r, ue, qw);

        d = u_h(6*(e-1) + (1:6)) - (M \ b);

        e2 = e2 + d.' * M * d;
        n2 = n2 + sum(sum(ue.^2, 1) .* qw);

    end

    eRel = sqrt(e2 / n2);

end

function [M, b] = rmGram(r, ue, qw)
% RM(E) mass matrix and load vector against a field, in the basis
% phi_i = e_i, phi_{3+i} = e_i x (x - x_E)

    nq = size(r,2);
    Phi = zeros(3, nq, 6);

    for i = 1:3
        ei = zeros(3,1); ei(i) = 1;
        Phi(i,:,i) = 1;
        Phi(:,:,3+i) = cross(repmat(ei,1,nq), r, 1);
    end

    M = zeros(6); b = zeros(6,1);

    for a = 1:6

        b(a) = sum(sum(Phi(:,:,a) .* ue, 1) .* qw);

        for c = a:6
            M(a,c) = sum(sum(Phi(:,:,a) .* Phi(:,:,c), 1) .* qw);
            M(c,a) = M(a,c);
        end

    end

end

function cd = patchCellData(cell_struct, face_struct, face_global_geom, ...
                            B_local, P_local, sigma_h, u_h, problem)
% per-cell fields for visualisation
%
%   displacement         cell mean of u_h. u_h|_E = alpha + omega x (x - x_E)
%                        and x_E is the centroid, so the mean is alpha
%   displacement_exact   cell mean of u
%   displacement_error   || Pi_RM u - u_h ||_{0,E}, the local Dassi E_u norm
%
%   stress               Pi_E sigma_h, constant on the element, row-major
%   stress_exact         cell mean of sigma
%   stress_error         || sigma - Pi_E sigma_h ||_{0,E}, the local Dassi
%                        E_sig,Pi norm
%
% the error fields are local L2 norms, so summing them in quadrature recovers
% the global Dassi norms: sqrt(sum(field.^2)) = || . ||_0 over the whole mesh

    nCells = numel(cell_struct);
    s = 1/sqrt(2);

    cd.displacement       = zeros(nCells, 3);
    cd.displacement_exact = zeros(nCells, 3);
    cd.displacement_error = zeros(nCells, 1);
    cd.stress             = zeros(nCells, 9);
    cd.stress_exact       = zeros(nCells, 9);
    cd.stress_error       = zeros(nCells, 1);

    for e = 1:nCells

        qp = cell_struct(e).quad_points;
        qw = cell_struct(e).quad_weights(:);
        Ve = sum(qw);
        r = qp - cell_struct(e).center(:);

        ue = problem.u(qp);

        cd.displacement(e,:) = u_h(6*(e-1) + (1:3)).';
        cd.displacement_exact(e,:) = (ue * qw / Ve).';

        % Pi_RM u - u_h, both expressed in the RM basis of this element
        d = u_h(6*(e-1) + (1:6)) - rmCoefficients(r, ue, qw);
        du = d(1:3) + cross(repmat(d(4:6),1,size(r,2)), r, 1);

        cd.displacement_error(e) = sqrt(sum(sum(du.^2, 1).' .* qw));

        sl = localStressDofs(e, cell_struct, face_struct, face_global_geom, B_local, sigma_h);

        c = P_local{e} * sl;
        Sh = [   c(1), s*c(6), s*c(5); ...
              s*c(6),    c(2), s*c(4); ...
              s*c(5), s*c(4),    c(3)];

        Sq = problem.sigma(qp);
        dS = Sq - Sh;

        cd.stress(e,:)       = reshape(Sh.', 1, 9);
        cd.stress_exact(e,:) = reshape(reshape(Sq, 9, []) * qw / Ve, 1, 9);
        cd.stress_error(e)   = sqrt(sum(reshape(sum(sum(dS.^2,1),2), [], 1) .* qw));

    end

end

function c = rmCoefficients(r, ue, qw)
% coefficients of the L2 projection of a field onto RM(E), in the basis
% phi_i = e_i, phi_{3+i} = e_i x (x - x_E)

    [M, b] = rmGram(r, ue, qw(:).');
    c = M \ b;

end

function s = ternary(c, a, b)
    if c, s = a; else, s = b; end
end

function setupPatchPath()

    root = fileparts(fileparts(mfilename('fullpath')));

    needed = {'geometry','operators','assembly','solvers','io','meshes', ...
              'verification_incompressibility'};

    for k = 1:numel(needed)
        d = fullfile(root, needed{k});
        if isfolder(d) && ~contains([path pathsep], [d pathsep])
            addpath(d);
        end
    end

end
