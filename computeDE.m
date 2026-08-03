function D_local = computeDE(cell_struct, A_local, B_local)
% construct the local divergence reconstruction matrix D_E

    nCells = numel(A_local);
    D_local = cell(nCells, 1);
    
    for e = 1:nCells
        A_E = A_local{e};
        % assert(rcond(A_E) > 1e-12, 'A_E is nearly singular.');

        volume = cell_struct(e).volume;
        B_E = B_local{e}.matrix;
    
        nCols = size(B_E, 2);
        D_E = zeros(6, nCols);
    
        for j = 1:nCols
    
            % translation part
            alpha = B_E(1:3, j) / volume;
    
            % rotational RHS
            b = B_E(4:6, j);
    
            % rotational coefficients
            omega = A_E \ b;
    
            D_E(:, j) = [alpha;omega];
        end
    
        D_local{e} = D_E;
    end

end