function [K_bc, RHS_bc, bc] = applyContactBoundaryConstraints(K_global, RHS, cell_struct, face_struct, face_global_geom)
% impose component-wise constant traction boundary conditions

    nTotal = size(K_global, 1);

    C = sparse(0, nTotal);
    d = zeros(0, 1);

    for f = 1:numel(face_struct)

        face = face_struct(f);

        if ~face.is_boundary
            continue;
        end

        mask = logical(face.traction_mask(:));

        if ~any(mask)
            continue;
        end

        geom = face_global_geom(f);

        % determine outward normal sign
        e = face.cells(face.cells > 0);
        outwardDirection = face.center(:) - cell_struct(e).center(:);
        normalSign = sign(dot(geom.n, outwardDirection));

        assert(normalSign ~= 0, 'Invalid boundary normal.');

        % construct constraints using the actual face basis
        qp = face.quad_points;
        w = face.quad_weights;

        M = zeros(6, 6);
        g = zeros(6, 1);

        P = diag(double(mask));
        traction = face.traction_value(:);

        for q = 1:numel(w)
            Phi = evaluateFaceBasis(qp(:,q), face.center(:), geom.Qf, geom.t1, geom.t2, geom.n, face.area, geom.m20, geom.m02);
            M = M + w(q) * (Phi.' * P * Phi);
            g = g + w(q) * (Phi.' * P * traction);
        end

        Cface = normalSign * M;

        % keep independent constraint equations
        [~, pivots] = rref(full(Cface));

        Cnew = sparse(numel(pivots), nTotal);
        Cnew(:, face.stress_dofs(:)) = Cface(pivots, :);

        C = [C; Cnew];
        d = [d; g(pivots)];

    end

    % enforce constraints using lagrange multipliers
    nConstraints = size(C, 1);

    K_bc = [K_global, C.'; C, sparse(nConstraints, nConstraints)];
    RHS_bc = [RHS; d];

    bc.nOriginal = nTotal;
    bc.nConstraints = nConstraints;

end