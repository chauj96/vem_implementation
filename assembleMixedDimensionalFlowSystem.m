function [A, blocks] = assembleMixedDimensionalFlowSystem( ...
    cell_struct, fracture, flowDofs, ...
    A_mfd_local, E_mfd_local, ...
    AF_mfd_local, EF_mfd_local, R, TF)

% Assemble the global mixed-dimensional Darcy system.
%
% Unknown ordering:
%   x = [w; p; wF; pF]
%
% Local MFD inner products use outward-oriented fluxes.
% Global flux DOFs use fixed face/edge orientations.
%
% Global matrix:
%
%   [ W+R    E'      0      TF' ]
%   [ -E     0       0       0  ]
%   [  0     0      WF      EF' ]
%   [ -TF    0      -EF      0  ]

    nW = flowDofs.nW;
    nP = flowDofs.nP;
    nWF = flowDofs.nWF;
    nPF = flowDofs.nPF;

    %% 1. Assemble bulk operators

    W = sparse(nW, nW);
    E = sparse(nP, nW);

    for e = 1:numel(cell_struct)

        dofs = cell_struct(e).flow_dofs(:);
        signs = cell_struct(e).flow_signs(:);

        S = diag(signs);

        Wlocal = S * A_mfd_local{e} * S;
        Elocal = E_mfd_local{e};

        pDof = cell_struct(e).pressure_dof;

        W(dofs, dofs) = W(dofs, dofs) + Wlocal;
        E(pDof, dofs) = E(pDof, dofs) + Elocal;
    end

    %% 2. Assemble fracture operators

    WF = sparse(nWF, nWF);
    EF = sparse(nPF, nWF);

    for k = 1:fracture.nFaces

        dofs = fracture.face_flux_dofs{k}(:);
        signs = fracture.face_flux_signs{k}(:);

        S = diag(signs);

        WFlocal = S * AF_mfd_local{k} * S;
        EFlocal = EF_mfd_local{k};

        pFDof = fracture.pressure_dofs(k);

        WF(dofs, dofs) = WF(dofs, dofs) + WFlocal;
        EF(pFDof, dofs) = EF(pFDof, dofs) + EFlocal;
    end

    %% 3. Add interface resistance

    Wstar = W + R;

    %% 4. Assemble global block matrix

    A = [Wstar, E', sparse(nW, nWF), TF';
         -E, sparse(nP, nP), sparse(nP, nWF), sparse(nP, nPF);
         sparse(nWF, nW), sparse(nWF, nP), WF, EF';
         -TF, sparse(nPF, nP), -EF, sparse(nPF, nPF)];

    %% 5. Store assembled blocks

    blocks.W = W;
    blocks.E = E;
    blocks.WF = WF;
    blocks.EF = EF;
    blocks.R = R;
    blocks.TF = TF;
    blocks.Wstar = Wstar;

    %% 6. Verification

    nTotal = nW + nP + nWF + nPF;

    assert(isequal(size(A), [nTotal, nTotal]), ...
        'Incorrect global matrix size.');

    symW = norm(W - W', 'fro') / max(norm(W, 'fro'), eps);
    symWF = norm(WF - WF', 'fro') / max(norm(WF, 'fro'), eps);
    symWstar = norm(Wstar - Wstar', 'fro') / max(norm(Wstar, 'fro'), eps);

    assert(symW < 1e-12, 'Bulk W is not symmetric.');
    assert(symWF < 1e-12, 'Fracture WF is not symmetric.');
    assert(symWstar < 1e-12, 'Wstar is not symmetric.');

    assert(nnz(E) == sum(arrayfun(@(c) numel(c.faces), cell_struct)), ...
        'Unexpected bulk divergence connectivity.');

    assert(nnz(EF) == sum(cellfun(@numel, fracture.face_edges)), ...
        'Unexpected fracture divergence connectivity.');

    assert(all(isfinite(nonzeros(A))), ...
        'Global matrix contains NaN or Inf.');

    fprintf('\n[# GLOBAL MIXED-DIMENSIONAL FLOW ASSEMBLY]\n');
    fprintf('  Bulk W: %d x %d\n', size(W, 1), size(W, 2));
    fprintf('  Bulk E: %d x %d\n', size(E, 1), size(E, 2));
    fprintf('  Fracture WF: %d x %d\n', size(WF, 1), size(WF, 2));
    fprintf('  Fracture EF: %d x %d\n', size(EF, 1), size(EF, 2));
    fprintf('  Global A: %d x %d\n', size(A, 1), size(A, 2));
    fprintf('  Global nonzeros: %d\n', nnz(A));
    fprintf('  Bulk W symmetry error: %.3e\n', symW);
    fprintf('  Fracture WF symmetry error: %.3e\n', symWF);
    fprintf('  Wstar symmetry error: %.3e\n', symWstar);
    fprintf('  Status: PASS\n');
end
