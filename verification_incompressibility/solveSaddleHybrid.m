function [sigma_h, u_h, info] = solveSaddleHybrid(cell_struct, face_struct, ...
        face_global_geom, B_local, K_cons, K_stab, rhs_g, rhs_f, tol, maxit)
% hybridized solve of the mixed saddle system
%
%   [ K  B'] [sigma]   [g]
%   [ B  0 ] [  u  ] = [f]
%
% The traction continuity built into Sigma_h is released: each element gets its
% own copy of the stress DOFs on its faces, and continuity is reimposed by a
% multiplier lambda on the interior faces. With K block diagonal, sigma_E and
% u_E are eliminated element by element and only lambda is solved globally:
%
%   [ K^  B^' C^'] [s^]   [g]
%   [ B^  0   0  ] [u ] = [f]      K^, B^ block diagonal
%   [ C^  0   0  ] [lam]   [0]
%
%   H lambda = h,   H = sum_E C_E S_E C_E',   S_E = (M_E^-1)_{ss}
%
% with M_E = [K_E B_E'; B_E 0] the local saddle matrix. H is symmetric positive
% (semi)definite of size 6*nInteriorFaces, so the indefinite MINRES solve is
% replaced by a CG solve on a smaller system.
%
% Because assembleGlobalMatrices maps every element into a common global face
% convention, the two copies of an interior face are compared directly and
% C_E is a signed selection: +1 on the first incident element, -1 on the second.
%
% sigma_h and u_h are returned in the same global convention and ordering as the
% monolithic path, so downstream code is unchanged.
%
% info.jump is the largest disagreement between the two copies of an interior
% face, the residual of the constraint the multiplier enforces

    if nargin < 9  || isempty(tol),   tol   = 1e-11; end
    if nargin < 10 || isempty(maxit), maxit = 2000;  end

    nCells = numel(cell_struct);
    nFaces = numel(face_struct);

    %% face -> cell incidence, and the multiplier numbering
    [faceCells, nAdj] = faceIncidence(cell_struct, nFaces);

    isInterior = (nAdj == 2);

    lamIdx = zeros(nFaces, 1);
    lamIdx(isInterior) = 1:nnz(isInterior);

    nL = 6 * nnz(isInterior);

    assert(nL > 0, 'No interior faces: nothing to hybridize.');

    %% pass 1: local elimination, assemble H and h
    %
    % T_E is diagonal +-1, so it is carried as a vector and reused in pass 2.
    % S_E is only needed on the interior-face rows and columns, so only those
    % columns of M_E^-1 are formed

    tCache = cell(nCells, 1);

    nH = 0;

    for e = 1:nCells
        nH = nH + (6 * nnz(lamIdx(cell_struct(e).faces)))^2;
    end

    IH = zeros(nH,1); JH = zeros(nH,1); VH = zeros(nH,1);
    h = zeros(nL,1);

    ptr = 0;

    t1 = tic;

    for e = 1:nCells

        [M, tE, sDofs, dDofs] = localSaddle(e, cell_struct, face_struct, ...
                                     face_global_geom, B_local, K_cons, K_stab);

        tCache{e} = tE;

        [loc, glb, sg] = interiorMaps(e, cell_struct, faceCells, lamIdx);

        n = numel(loc);
        nl = size(M,1);

        % one factorization per element: the lambda = 0 solve and the columns of
        % M_E^-1 on the interior-face DOFs are solved against the same matrix
        E = zeros(nl, n+1);
        E(:,1) = [rhs_g(sDofs); rhs_f(dDofs)];
        E(sub2ind([nl n+1], loc(:), (2:n+1).')) = 1;

        W = solveLocalSaddle(M, E);

        if n == 0
            continue;
        end

        h(glb) = h(glb) + sg .* W(loc,1);

        Sloc = (sg .* W(loc,2:end)) .* sg.';   % C_E S_E C_E'

        [rows, cols] = ndgrid(glb, glb);

        idx = ptr + (1:n^2);

        IH(idx) = rows(:);
        JH(idx) = cols(:);
        VH(idx) = Sloc(:);

        ptr = ptr + n^2;

    end

    tPass1 = toc(t1);

    t2 = tic;

    H = sparse(IH(1:ptr), JH(1:ptr), VH(1:ptr), nL, nL);

    clear IH JH VH

    % H is assembled from symmetric blocks; enforce exact symmetry for ichol
    H = (H + H.') / 2;

    %% solve the multiplier system
    L = robustIchol(H);

    [lambda, flag, relres, iter] = pcg(H, h, tol, maxit, L, L.');

    tSolve = toc(t2);

    t3 = tic;

    %% pass 2: recover sigma and u element by element
    sigma_h = zeros(6*nFaces, 1);
    u_h = zeros(6*nCells, 1);

    seen = false(6*nFaces, 1);
    jump = 0;

    for e = 1:nCells

        [M, ~, sDofs, dDofs] = localSaddle(e, cell_struct, face_struct, ...
                                   face_global_geom, B_local, K_cons, K_stab, tCache{e});

        [loc, glb, sg] = interiorMaps(e, cell_struct, faceCells, lamIdx);

        rhsE = [rhs_g(sDofs); rhs_f(dDofs)];
        rhsE(loc) = rhsE(loc) - sg .* lambda(glb);

        w = solveLocalSaddle(M, rhsE);

        sE = w(1:numel(sDofs));

        % the two copies of an interior face agree only up to the CG tolerance;
        % compare against whatever the first element wrote, then overwrite
        old = seen(sDofs);

        if any(old)
            jump = max(jump, max(abs(sigma_h(sDofs(old)) - sE(old))));
        end

        sigma_h(sDofs) = sE;
        seen(sDofs) = true;

        u_h(dDofs) = w(numel(sDofs)+1:end);

    end

    info.flag = flag;
    info.iterations = iter;
    info.relres = relres;
    info.multipliers = nL;
    info.jump = jump;

    info.tPass1 = tPass1;      % local elimination and assembly of H
    info.tSolve = tSolve;      % ichol and pcg on the multiplier system
    info.tPass2 = toc(t3);     % elementwise recovery of sigma and u

end

%% ------------------------------------------------------------------

function X = solveLocalSaddle(M, RHS)
% solve M X = RHS after symmetric equilibration, returning X in the original
% variables
%
%   M x = b   <=>   (S M S)(S^-1 x) = S b,   S = diag(s)
%
% The zero (2,2) block makes the raw local saddle matrix badly scaled: on the
% two-fault mesh cond(M_E) reaches 2.4e25 and MATLAB reports rcond ~ 1e-19,
% while row-inf-norm equilibration brings the worst case to 4.4e8. B_E has full
% row rank on every element, so this is scaling and not a degeneracy.

    s = 1 ./ sqrt(max(max(abs(M), [], 2), realmin));

    % (s*s') .* M is bitwise symmetric, so the solver keeps the LDL path;
    % (s.*M).*s' is not, and falls back to general LU
    X = s .* (((s * s.') .* M) \ (s .* RHS));

end

function [M, tE, sDofs, dDofs] = localSaddle(e, cell_struct, face_struct, ...
                                     face_global_geom, B_local, K_cons, K_stab, tE)
% local saddle matrix in the global face convention,
%   M_E = [K_E B_E'; B_E 0],  K_E = T_E' (Kc + Ks) T_E,  B_E = B T_E

    fl = cell_struct(e).faces(:).';
    nf = numel(fl);

    sDofs = reshape(6*(fl-1) + (1:6).', 1, []);
    dDofs = 6*(e-1) + (1:6);

    if nargin < 8 || isempty(tE)

        tE = zeros(6*nf, 1);

        for lf = 1:nf

            f = fl(lf);

            tE(6*(lf-1) + (1:6)) = faceSignTransformation(f, B_local{e}.geom(lf), ...
                                                          face_global_geom(f));

        end

    end

    KE = (tE .* (K_cons{e} + K_stab{e})) .* tE.';
    BE = B_local{e}.matrix .* tE.';

    M = [KE, BE.'; BE, zeros(6)];

end

function [loc, glb, sg] = interiorMaps(e, cell_struct, faceCells, lamIdx)
% local stress DOFs on interior faces, the multiplier DOFs they constrain, and
% the sign of C_E: +1 on the first incident element, -1 on the second

    fl = cell_struct(e).faces(:).';

    keep = find(lamIdx(fl).' > 0);

    loc = reshape(6*(keep-1) + (1:6).', [], 1);
    glb = reshape(6*(lamIdx(fl(keep)).'-1) + (1:6).', [], 1);

    s = ones(1, numel(keep));
    s(faceCells(fl(keep), 1).' ~= e) = -1;

    sg = reshape(repmat(s, 6, 1), [], 1);

end

function [faceCells, nAdj] = faceIncidence(cell_struct, nFaces)
% face -> cell incidence rebuilt from the cell face lists

    faceCells = zeros(nFaces, 2);
    nAdj = zeros(nFaces, 1);

    for e = 1:numel(cell_struct)

        for f = cell_struct(e).faces(:).'

            nAdj(f) = nAdj(f) + 1;

            assert(nAdj(f) <= 2, 'Face %d has more than two incident cells.', f);

            faceCells(f, nAdj(f)) = e;

        end

    end

end

function L = robustIchol(M)
% incomplete Cholesky, retried with a growing diagonal shift on breakdown

    shift = 0;

    for attempt = 1:8
        try
            L = ichol(M, struct('type','ict','droptol',1e-3,'diagcomp',shift));
            return;
        catch
            if shift == 0
                shift = 1e-4;
            else
                shift = shift * 10;
            end
        end
    end

    L = spdiags(sqrt(max(full(diag(M)), eps)), 0, size(M,1), size(M,2));

end
