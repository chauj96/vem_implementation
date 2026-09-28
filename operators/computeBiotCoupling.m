function D_biot_local = computeBiotCoupling(cell_struct, P_local)
% construct the local Biot stress-pressure coupling matrix D_E
%
% D_E = alpha * |E| * (A I)^T * P_E

    nCells = numel(cell_struct);
    D_biot_local = cell(nCells, 1);

    % identity tensor in the symmetric tensor basis
    % ordering: xx, yy, zz, yz, xz, xy
    Ivec = [1; 1; 1; 0; 0; 0];

    for e = 1:nCells

        volume = cell_struct(e).volume;
        alpha = cell_struct(e).alpha;
        A = cell_struct(e).Cinv;

        P_E = P_local{e};

        % coefficient vector of A I
        aI = A * Ivec;

        % local Biot coupling: 1 x nStressDofs_E
        D_biot_local{e} = alpha * volume * (aI.' * P_E);

    end

end