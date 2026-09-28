function E_mfd_local = computeMFDCoupling(cell_struct)

    nCells = numel(cell_struct);
    E_mfd_local = cell(nCells, 1);

    for e = 1:nCells
        s = cell_struct(e).faces_orientation(:).';
        E_mfd_local{e} = -s;
    end

end