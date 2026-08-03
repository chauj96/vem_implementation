function A_local = computeAEQuadrature(cell_struct)
% construct the inertia matrix A_E by quadrature rule

    nCells = numel(cell_struct);
    A_local = cell(nCells, 1);
    
    for e = 1:nCells
        xE = cell_struct(e).center(:);  % barycenter of current cell
    
        qp = cell_struct(e).quad_points;
        w = cell_struct(e).quad_weights;
        nQuad = cell_struct(e).nQuad;
    
        A_E = zeros(3,3);
    
        for q = 1:nQuad
            r = qp(:, q) - xE;
            A_E = A_E + w(q) * ((r' * r) * eye(3) - r * r');
        end
    
        A_local{e} = A_E;
    end

end