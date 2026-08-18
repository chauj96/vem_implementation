function rhs_f = assembleBodyForce(cell_struct, problem)
% assemble the body force contribution to the displacement RHS
%
%   (rhs_f)_{E,j} = integral_E f . R_j dV
%
% with the rigid-motion basis of cell E
%
%   R_1..R_3 = e_1..e_3,        R_{3+k} = e_k x (x - x_E)
%
% so that the rotational entries are  e_k . integral_E (x - x_E) x f dV.
% this is the right-hand side of the discrete divergence equation B sigma = rhs_f,
% matching div(sigma) = f

    nCells = numel(cell_struct);
    rhs_f = zeros(6*nCells, 1);

    for e = 1:nCells

        qp = cell_struct(e).quad_points;
        qw = cell_struct(e).quad_weights(:).';

        xE = cell_struct(e).center(:);

        fq = problem.f(qp);                 % 3 x nq

        % translations
        trans = sum(fq .* qw, 2);

        % rotations: integral (x - xE) x f
        r = qp - xE;
        rot = sum(cross(r, fq, 1) .* qw, 2);

        rhs_f(6*(e-1) + (1:6)) = [trans; rot];

    end

end
