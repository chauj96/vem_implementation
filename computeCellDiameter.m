function cell_struct = computeCellDiameter(cell_struct, V3, cells3D)
% precompute the diameter of every polyhedral element

    nCells = numel(cell_struct);

    for e = 1:nCells
        
        % going to vectorize the computation
        X = V3(cells3D{e}, :); % collect whole vertices of current cell
        D = pdist(X);

        cell_struct(e).diameter = max(D);

    end
end