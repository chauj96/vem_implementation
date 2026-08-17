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
        
        % precompute cell diameter (kinda vectorization)
        X = V3(cells3D{e}, :); 
        pairwiseDist = pdist(X);

        cell_struct(e).diameter = max(pairwiseDist);

        % assign the compliance tensor
        cell_struct(e).Cinv = Cinv;

    end

end