function [sigma_h, u_h] = solveHybridVEM(cell_struct, face_struct, face_global_geom, B_local, K_cons, K_stab, rhs_D, rhs_f, lambdaMat, muMat, tol, maxit)
% solve the hybridized mixed-VEM system
%
% element stress DOFs are kept independent and traction continuity on
% interior faces is enforced through the multiplier lambda
%
% Local system:
%
%   M_E [sigma_E; u_E] - C_E' lambda = b_E
%
% where
%
%   M_E = [K_E  B_E';
%          B_E   0   ].
%
% element unknowns are eliminated locally giving the global multiplier system
%
%   H lambda = h.
%
% once lambda is obtained, sigma_E and u_E are recovered independently on each element

    if nargin < 11 || isempty(tol)
        tol = 1e-11;
    end

    if nargin < 12 || isempty(maxit)
        maxit = 2000;
    end

    nCells = numel(cell_struct);
    nFaces = numel(face_struct);

    % face incidence and multiplier numbering
    [faceCells, nAdj] = faceIncidence(cell_struct, nFaces);

    isInterior = (nAdj == 2);

    lamIdx = zeros(nFaces,1);
    lamIdx(isInterior) = 1:nnz(isInterior);

    nL = 6 * nnz(isInterior);

    assert(nL > 0, 'No interior faces: nothing to hybridize.');


    % compute and assign neumann data
    [neuDofs, neuVals] = computeNeumannData(face_struct, face_global_geom, lambdaMat, muMat);

    neuMask = false(6*nFaces,1);
    neuData = zeros(6*nFaces,1);

    neuMask(neuDofs) = true;
    neuData(neuDofs) = neuVals;


    %%  PASS 1: local elimination and assembly of H lambda = h
    tCache = cell(nCells,1);

    % estimate number of entries required for H
    nH = 0;

    for e = 1:nCells

        elementFaces = cell_struct(e).faces;

        nLocalLambdaDofs = 6 * nnz(lamIdx(elementFaces));

        nH = nH + nLocalLambdaDofs^2;

    end

    IH = zeros(nH,1);
    JH = zeros(nH,1);
    VH = zeros(nH,1);

    h = zeros(nL,1);

    ptr = 0;

    t1 = tic;

    for e = 1:nCells

        % local saddle matrix in global face convention
        [M, TE, sDofs, dDofs] = localSaddle(e, cell_struct, face_struct, face_global_geom, B_local, K_cons, K_stab);

        tCache{e} = TE;

        % local/global multiplier maps
        [loc, glb, sg] = interiorMaps(e, cell_struct, faceCells, lamIdx);

        % element RHS
        rhsE = [rhs_D(sDofs);
                rhs_f(dDofs)];

        % impose prescribed Neumann stress DOFs locally
        [M, rhsE] = applyLocalNeumann(M, rhsE, sDofs, neuMask, neuData);

        n = numel(loc);
        nl = size(M,1);

        % first column: M_E^{-1} b_E
        % remaining columns: response to the interior-face multiplier DOFs
        R = zeros(nl,n+1);

        R(:,1) = rhsE;

        if n > 0
            R(sub2ind([nl,n+1], loc(:), (2:n+1).')) = 1;
        end

        W = solveLocalSaddle(M,R);

        if n == 0
            continue;
        end

        % contribution to h
        h(glb) = h(glb) + sg .* W(loc,1);

        % contribution to H
        Sloc = (sg .* W(loc,2:end)) .* sg.';

        [rows,cols] = ndgrid(glb,glb);

        idx = ptr + (1:n^2);

        IH(idx) = rows(:);
        JH(idx) = cols(:);
        VH(idx) = Sloc(:);

        ptr = ptr + n^2;

    end

    tPass1 = toc(t1);

    %  assemble and solve multiplier system
    fprintf('\n[# SOLVE THE LINEAR SYSTEM OF HYBRIDIZED VEM]\n');

    H = sparse(IH(1:ptr), JH(1:ptr), VH(1:ptr), nL, nL);

    clear IH JH VH

    % remove small numerical asymmetry
    H = (H + H.') / 2;

    t2 = tic;

    % incomplete Cholesky preconditioner
    L = robustIchol(H);

    [lambda_h, flag, relres, iter] = pcg(H, h, tol, maxit, L, L.');

    tSolve = toc(t2);

    %%  PASS 2: elementwise recovery of sigma (stress) AND u (displacement)
    t3 = tic;

    sigma_h = zeros(6*nFaces,1);
    u_h = zeros(6*nCells,1);

    seen = false(6*nFaces,1);

    jump = 0;

    for e = 1:nCells

        [M, ~, sDofs, dDofs] = localSaddle(e, cell_struct, face_struct, face_global_geom, B_local, K_cons, K_stab, tCache{e});

        [loc, glb, sg] = interiorMaps(e, cell_struct, faceCells, lamIdx);

        rhsE = [rhs_D(sDofs);
                rhs_f(dDofs)];

        % impose Neumann stress DOFs
        [M, rhsE] = applyLocalNeumann(M, rhsE, sDofs, neuMask, neuData);

        % multiplier contribution
        rhsE(loc) = rhsE(loc) - sg .* lambda_h(glb);

        % recover local stress and displacement
        w = solveLocalSaddle(M,rhsE);

        nLocalStressDofs = numel(sDofs);

        sigmaE = w(1:nLocalStressDofs);
        uE = w(nLocalStressDofs+1:end);

        % compare the two copies on interior faces
        old = seen(sDofs);

        if any(old)
            jump = max(jump, max(abs(sigma_h(sDofs(old)) - sigmaE(old))));
        end

        % return solution in standard global ordering
        sigma_h(sDofs) = sigmaE;
        seen(sDofs) = true;

        u_h(dDofs) = uE;

    end

    tPass2 = toc(t3);

    %  SOLVER INFORMATION
    % info.flag = flag;
    % info.iterations = iter;
    % info.relres = relres;
    % info.multipliers = nL;
    % info.jump = jump;
    % 
    % info.tPass1 = tPass1;
    % info.tSolve = tSolve;
    % info.tPass2 = tPass2;
end

%% helper functions

% local saddle matrix
function [M, TE, sDofs, dDofs] = localSaddle(e, cell_struct, face_struct, face_global_geom, B_local, K_cons, K_stab, TE)

    elementFaces = cell_struct(e).faces(:).';

    nElementFaces = numel(elementFaces);
    nLocalStressDofs = 6 * nElementFaces;

    % global face DOF indices used only for accessing global-face data
    sDofs = reshape(6*(elementFaces-1) + (1:6).', 1, []);

    dDofs = 6 * (e-1) + (1:6);

    % transformation from local face convention to common global face convention
    if nargin < 8 || isempty(TE)

        TE = zeros(nLocalStressDofs);

        for lf = 1:nElementFaces

            f = elementFaces(lf);

            localIdx = 6 * (lf-1) + (1:6);

            localGeom = B_local{e}.geom(lf);
            globalGeom = face_global_geom(f);

            Tf = computeFaceTransformation(f, localGeom, globalGeom, face_struct);

            assert(norm(Tf' * Tf - eye(6),'fro') < 1e-12, 'Invalid face transformation: cell %d, face %d.', e, f);

            TE(localIdx,localIdx) = Tf;

        end

    end

    % transform local VEM matrices to global face convention
    KE = TE' * (K_cons{e} + K_stab{e}) * TE;

    BE = B_local{e}.matrix * TE;

    % local saddle matrix
    M = [KE, BE';
         BE, zeros(6)];

end


% apply neumann condition locally
function [M, rhsE] = applyLocalNeumann(M, rhsE, sDofs, neuMask, neuData)

    isNeu = neuMask(sDofs);

    if ~any(isNeu)
        return;
    end

    % positions inside the element-local stress vector
    localNeu = find(isNeu);

    % prescribed traction DOFs
    neuVals = neuData(sDofs(localNeu));

    % move known stress contributions to RHS
    rhsE = rhsE - M(:,localNeu) * neuVals;

    % impose sigma_N = prescribed value
    M(:,localNeu) = 0;
    M(localNeu,:) = 0;

    M(localNeu,localNeu) = eye(numel(localNeu));

    rhsE(localNeu) = neuVals;

end


% interior face multiplier map 
function [loc, glb, sg] = interiorMaps(e, cell_struct, faceCells, lamIdx)

    elementFaces = cell_struct(e).faces(:).';

    % interior faces belonging to element e
    keep = find(lamIdx(elementFaces).' > 0);

    % positions in the element-local stress vector
    loc = reshape(6*(keep-1) + (1:6).', [], 1);

    % corresponding global multiplier DOFs
    glb = reshape(6*(lamIdx(elementFaces(keep)).'-1) + (1:6).', [], 1);

    % C_E orientation:
    % +I for first incident element
    % -I for second incident element
    s = ones(1,numel(keep));

    s(faceCells(elementFaces(keep),1).' ~= e) = -1;

    sg = reshape(repmat(s,6,1), [],1);

end


% face -> cell incidence
function [faceCells, nAdj] = faceIncidence(cell_struct, nFaces)

    faceCells = zeros(nFaces,2);
    nAdj = zeros(nFaces,1);

    for e = 1:numel(cell_struct)

        for f = cell_struct(e).faces(:).'

            nAdj(f) = nAdj(f) + 1;

            assert(nAdj(f) <= 2, 'Face %d has more than two incident cells.', f);

            faceCells(f,nAdj(f)) = e;

        end

    end

end


% local saddle solve
function X = solveLocalSaddle(M,RHS)

    % symmetric row/column equilibration
    s = 1 ./ sqrt(max(max(abs(M),[],2),realmin));

    X = s .* (((s*s.') .* M) \ (s .* RHS));

end


% robust incomplete cholesky
function L = robustIchol(H)

    shift = 0;

    for attempt = 1:8

        try
            L = ichol(H, struct('type','ict', 'droptol',1e-3, 'diagcomp',shift));
            return;

        catch

            if shift == 0
                shift = 1e-4;
            else
                shift = 10*shift;
            end

        end

    end

    % diagonal fallback
    L = spdiags(sqrt(max(full(diag(H)),eps)), 0, size(H,1), size(H,2));

end