function K_cons = computeConsistency(cell_struct, P_local)
% construct the local consistency matrix K_E^{cons}

    nCells = numel(cell_struct);
    K_cons = cell(nCells,1);

    for e = 1:nCells

        volume = cell_struct(e).volume;
        Cinv = cell_struct(e).Cinv;
        P = P_local{e};

        K_cons{e} = volume * (P.' * Cinv * P);

    end

end