function P_local = computeProjectionMatrix(cell_struct, Bproj)
% Construct the polynomial projection matrix P_E for every element.
% M_E P_E = Bproj_E where M_E = |E| I_6.

    nCells = numel(cell_struct);
    P_local = cell(nCells, 1);

    for e = 1:nCells

        volume = cell_struct(e).volume;
        Bproj_E = Bproj{e}.matrix;

        P_local{e} = Bproj_E / volume;
    end
end