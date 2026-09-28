function M_biot_local = computeBiotStorage(cell_struct)
% construct the local Biot pressure storage matrix M_E
%
% M_E = |E| * (s0 + alpha^2 / kappa)

    nCells = numel(cell_struct);
    M_biot_local = cell(nCells, 1);

    for e = 1:nCells

        volume = cell_struct(e).volume;
        s0     = cell_struct(e).s0;
        alpha  = cell_struct(e).alpha;
        kappa  = cell_struct(e).kappa;

        M_biot_local{e} = volume * (s0 + alpha^2 / kappa);

    end

end