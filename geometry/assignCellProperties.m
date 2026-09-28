function cell_struct = assignCellProperties(cell_struct, V3, cells3D, lambda, mu, alpha, s0)
% assign cell-level geometric and material properties
% (1) element diameter, (2) compliance tensor

    % compute compliance tensor C^{-1} (for now assume isotropic)
    a = (lambda + mu) / (mu * (3*lambda + 2*mu));
    b = -lambda / (2*mu * (3*lambda + 2*mu));
    c = 1 / (2*mu);
    kappa = lambda + 2*mu/3;
    
    Cinv = [ ...
        a, b, b, 0, 0, 0;
        b, a, b, 0, 0, 0;
        b, b, a, 0, 0, 0;
        0, 0, 0, c, 0, 0;
        0, 0, 0, 0, c, 0;
        0, 0, 0, 0, 0, c];

    % flow properties
    K_perm = eye(3);

    nCells = numel(cell_struct);

    for e = 1:nCells
        
        % precompute cell diameter (kinda vectorization)
        X = V3(cells3D{e}, :);

        % pairwise squared distances via the Gram matrix
        % (avoids pdist, which needs the Statistics Toolbox)
        sq = sum(X.^2, 2);
        D2 = sq + sq.' - 2*(X*X.');

        cell_struct(e).diameter = sqrt(max(max(D2(:)), 0));

        % assign the compliance tensor
        cell_struct(e).Cinv = Cinv;

        cell_struct(e).alpha = alpha;
        cell_struct(e).kappa = kappa;
        cell_struct(e).s0 = s0;
        cell_struct(e).lambda = lambda;
        cell_struct(e).mu = mu;

        % flow
        cell_struct(e).perm = K_perm;

    end

end