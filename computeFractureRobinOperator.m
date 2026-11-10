function R = computeFractureRobinOperator(fracture, face_struct, flowDofs, Tminus, Tplus)
% construct the bulk-fracture Robin resistance matrix.
%
% interface law:
%   lambda^pm = T^pm * (p^pm - pF)
%
% integrated normal flux DOFs:
%   q_f = integral_f (w . n) dS
%
% R contributes to the bulk flux block:
%   Wstar = W + R
%
% Tminus and Tplus are transmissibilities per unit area.
% They may be scalars or vectors of length fracture.nFaces.

    nF = fracture.nFaces;
    nW = flowDofs.nW;

    if isscalar(Tminus)
        Tminus = Tminus * ones(nF, 1);
    end

    if isscalar(Tplus)
        Tplus = Tplus * ones(nF, 1);
    end

    Tminus = Tminus(:);
    Tplus = Tplus(:);

    assert(numel(Tminus) == nF, 'Invalid minus-side transmissibility size.');
    assert(numel(Tplus) == nF, 'Invalid plus-side transmissibility size.');
    assert(all(isfinite(Tminus)) && all(Tminus > 0), 'Invalid Tminus.');
    assert(all(isfinite(Tplus)) && all(Tplus > 0), 'Invalid Tplus.');

    R = sparse(nW, nW);

    for k = 1:nF

        f = fracture.faces(k);
        Af = face_struct(f).area;

        dofMinus = face_struct(f).flow_dof;
        dofPlus = face_struct(f).flow_dof_plus;

        R(dofMinus, dofMinus) = R(dofMinus, dofMinus) + 1 / (Tminus(k) * Af);
        R(dofPlus, dofPlus) = R(dofPlus, dofPlus) + 1 / (Tplus(k) * Af);
    end

    % %% Verification
    % 
    % sym_err = norm(R - R.', 'fro') / max(norm(R, 'fro'), eps);
    % 
    % expectedNnz = 2 * nF;
    % 
    % assert(nnz(R) == expectedNnz, 'Unexpected number of Robin entries.');
    % assert(sym_err < 1e-12, 'Robin matrix is not symmetric.');
    % assert(all(nonzeros(R) > 0), 'Robin resistance must be positive.');
    % 
    % fprintf('\n[# BULK-FRACTURE ROBIN OPERATOR]\n');
    % fprintf('  Matrix size: %d x %d\n', size(R, 1), size(R, 2));
    % fprintf('  Nonzero entries: %d\n', nnz(R));
    % fprintf('  Symmetry error: %.3e\n', sym_err);
    % fprintf('  Status: PASS\n');
end
