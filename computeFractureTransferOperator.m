function TF = computeFractureTransferOperator(fracture, cell_struct, face_struct, flowDofs)
% construct the bulk-fracture transfer operator.
%
% TF maps globally oriented bulk flux DOFs to integrated flux entering the fracture.
%
% (TF * w)_k = q_minus(k) + q_plus(k)
%
% q_minus and q_plus use outward normals of the
% corresponding bulk cells.
%
% TF size: nPF x nW

    nF = fracture.nFaces;
    nW = flowDofs.nW;
    nPF = flowDofs.nPF;

    TF = sparse(nPF, nW);

    for k = 1:nF

        f = fracture.faces(k);

        em = fracture.minus_cells(k);
        ep = fracture.plus_cells(k);

        lfm = find(cell_struct(em).faces == f);
        lfp = find(cell_struct(ep).faces == f);

        assert(numel(lfm) == 1, 'Invalid minus-side face connectivity.');
        assert(numel(lfp) == 1, 'Invalid plus-side face connectivity.');

        sm = cell_struct(em).faces_orientation(lfm);
        sp = cell_struct(ep).faces_orientation(lfp);

        assert(sm == -sp, 'Inconsistent fracture normal orientation.');

        dofMinus = face_struct(f).flow_dof;
        dofPlus = face_struct(f).flow_dof_plus;

        pFDof = fracture.pressure_dofs(k);

        TF(pFDof, dofMinus) = TF(pFDof, dofMinus) + sm;
        TF(pFDof, dofPlus) = TF(pFDof, dofPlus) + sp;
    end

    %% Verification

    assert(all(sum(spones(TF), 2) == 2), 'Each fracture face must couple to two bulk flux DOFs.');

    fprintf('\n[# BULK-FRACTURE TRANSFER OPERATOR]\n');
    fprintf('  Matrix size: %d x %d\n', size(TF, 1), size(TF, 2));
    fprintf('  Nonzero entries: %d\n', nnz(TF));
    fprintf('  Status: PASS\n');
end
