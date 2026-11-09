function test_flowSignConvention(cell_struct, fracture, E_mfd_local, EF_mfd_local, TF)

    nCells = numel(cell_struct);
    nF = fracture.nFaces;

    %% Check bulk divergence convention

    for e = 1:nCells

        s = cell_struct(e).faces_orientation(:)';

        E = E_mfd_local{e};

        assert(isequal(E, -s), 'Bulk divergence sign mismatch in cell %d.', e);
    end

    %% Check fracture divergence convention

    for k = 1:nF

        s = fracture.face_edge_signs{k}(:)';

        EF = EF_mfd_local{k};

        assert(isequal(EF, -s), 'Fracture divergence sign mismatch in face %d.', k);
    end

    %% Check interface transfer convention

    for k = 1:nF

        f = fracture.faces(k);

        em = fracture.minus_cells(k);
        ep = fracture.plus_cells(k);

        lfm = find(cell_struct(em).faces == f);
        lfp = find(cell_struct(ep).faces == f);

        sm = cell_struct(em).faces_orientation(lfm);
        sp = cell_struct(ep).faces_orientation(lfp);

        dofMinus = cell_struct(em).flow_dofs(lfm);
        dofPlus = cell_struct(ep).flow_dofs(lfp);

        row = full(TF(fracture.pressure_dofs(k), :));

        assert(row(dofMinus) == sm, 'Minus-side transfer sign mismatch.');

        assert(row(dofPlus) == sp, 'Plus-side transfer sign mismatch.');

        assert(nnz(row) == 2, 'Unexpected transfer connectivity.');
    end

    fprintf('\n[# FLOW SIGN CONVENTION VERIFICATION]\n');
    fprintf('  Bulk divergence signs: PASS\n');
    fprintf('  Fracture divergence signs: PASS\n');
    fprintf('  Interface transfer signs: PASS\n');
    fprintf('  Status: PASS\n');
end
