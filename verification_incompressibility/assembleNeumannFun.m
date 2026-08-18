function [A_global, rhs, neuDofs, neuVals] = assembleNeumannFun(face_struct, A_global, rhs, face_global_geom, problem)
% impose prescribed tractions taken from an analytical solution
%
% variant of assembleNeumann for a general (non-constant) exact stress. on each
% Neumann face the exact traction is projected onto the six face basis
% functions by solving the 6x6 face mass system
%
%   M_ij = int_f phi_i . phi_j ,   b_i = int_f phi_i . (sigma n) ,   c = M \ b
%
% the resulting DOFs are then eliminated symmetrically, as in assembleNeumann

    neumannFaces = find([face_struct.is_neumann]);

    nNeu = numel(neumannFaces);
    neuDofs = zeros(1, 6*nNeu);
    neuVals = zeros(6*nNeu, 1);

    for k = 1:nNeu

        f = neumannFaces(k);

        xf = face_struct(f).center(:);
        Af = face_struct(f).area;
        qp = face_struct(f).quad_points;
        w = face_struct(f).quad_weights(:);

        g = face_global_geom(f);

        % exact traction with respect to the reference normal
        S = problem.sigma(qp);
        T = reshape(sum(S .* reshape(g.n, 1, 3, []), 2), 3, []);

        M = zeros(6,6);
        b = zeros(6,1);

        for q = 1:numel(w)
            Phi = evaluateFaceBasis(qp(:,q), xf, g.Qf, g.t1, g.t2, g.n, Af, g.m20, g.m02);
            M = M + w(q) * (Phi.' * Phi);
            b = b + w(q) * (Phi.' * T(:,q));
        end

        neuDofs(6*(k-1) + (1:6)) = 6*(f-1) + (1:6);
        neuVals(6*(k-1) + (1:6)) = M \ b;

    end

    if isempty(neuDofs)
        return;
    end

    % move the known contributions to the right-hand side, then eliminate
    rhs = rhs - A_global(:, neuDofs) * neuVals;

    A_global(:, neuDofs) = 0;
    A_global(neuDofs, :) = 0;
    A_global(neuDofs, neuDofs) = speye(numel(neuDofs));

    rhs(neuDofs) = neuVals;

end
