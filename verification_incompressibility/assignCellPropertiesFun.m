function cell_struct = assignCellPropertiesFun(cell_struct, V3, cells3D, problem)
% assign cell diameter and a heterogeneous compliance tensor
%
% variant of assignCellProperties for a spatially varying medium: mu and
% lambda are sampled from problem at the cell barycenter, so the material is
% piecewise constant per cell as in section 5.1.1

    nCells = numel(cell_struct);

    for e = 1:nCells

        % h_E = max_{i,j} |x_i - x_j|,  with
        % |x_i - x_j|^2 = |x_i|^2 + |x_j|^2 - 2 x_i . x_j
        X = V3(cells3D{e}, :);

        sq = sum(X.^2, 2);
        pairwiseDist2 = sq + sq.' - 2*(X*X.');

        cell_struct(e).diameter = sqrt(max(max(pairwiseDist2, 0), [], 'all'));

        xE = cell_struct(e).center(:);

        mu = problem.mu(xE);
        lambda = problem.lambda(xE);

        a = (lambda + mu) / (mu * (3*lambda + 2*mu));
        b = -lambda / (2*mu * (3*lambda + 2*mu));
        c = 1 / (2*mu);

        cell_struct(e).Cinv = [ ...
            a, b, b, 0, 0, 0;
            b, a, b, 0, 0, 0;
            b, b, a, 0, 0, 0;
            0, 0, 0, c, 0, 0;
            0, 0, 0, 0, c, 0;
            0, 0, 0, 0, 0, c];

        cell_struct(e).mu = mu;
        cell_struct(e).lambda = lambda;

    end

end
