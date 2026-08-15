function cell_struct = assignCellProperties(cell_struct, V3, cells3D, lambda, mu)
% assign cell-level geometric and material properties
% (1) element diameter, (2) compliance tensor

    % compute compliance tensor C^{-1} (for now assume isotropic)
    a = (lambda + mu) / (mu * (3*lambda + 2*mu));
    b = -lambda / (2*mu * (3*lambda + 2*mu));
    c = 1 / (2*mu);
    
    Cinv = [ ...
        a, b, b, 0, 0, 0;
        b, a, b, 0, 0, 0;
        b, b, a, 0, 0, 0;
        0, 0, 0, c, 0, 0;
        0, 0, 0, 0, c, 0;
        0, 0, 0, 0, 0, c];

    nCells = numel(cell_struct);

    for e = 1:nCells
        
        % h_E = max_{i,j} |x_i - x_j|,  with
        % |x_i - x_j|^2 = |x_i|^2 + |x_j|^2 - 2 x_i . x_j
        X = V3(cells3D{e}, :);

        sq = sum(X.^2, 2);
        pairwiseDist2 = sq + sq.' - 2*(X*X.');

        cell_struct(e).diameter = sqrt(max(max(pairwiseDist2, 0), [], 'all'));

        % assign the compliance tensor
        cell_struct(e).Cinv = Cinv;

    end

end