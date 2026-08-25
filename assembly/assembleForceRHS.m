function rhs_f = assembleForceRHS(cell_struct, bodyForce)
% assemble body-force contribution to the displacement RHS
%
%   (rhs_f)_{E,j} = integral_E f . R_j dV
%
% where R_j are the six rigid-motion basis functions on element E
%
% NOTE:
% This function computes the projection of f onto RM(E).
% For the convention -div(sigma) = f, use -rhs_f in the
% global saddle-point RHS.

    nCells = numel(cell_struct);
    rhs_f = zeros(6*nCells, 1);

    for e = 1:nCells

        qp = cell_struct(e).quad_points;
        qw = cell_struct(e).quad_weights(:).';

        xE = cell_struct(e).center(:);

        % evaluate body force at quadrature points
        fq = bodyForce(qp);

        % translational rigid-motion modes
        trans = sum(fq .* qw, 2);

        % rotational rigid-motion modes
        r = qp - xE;
        rot = sum(cross(r, fq, 1) .* qw, 2);

        rhs_f(6*(e-1) + (1:6)) = [trans; rot];

    end

end