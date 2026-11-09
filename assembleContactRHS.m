function RHS = assembleContactRHS(cell_struct, face_struct, face_global_geom, nStressDofs, nDispDofs, nContactTractionDofs, nInterfaceDispDofs)
% assemble RHS for the inclined fracture compression problem.
%
% dirichlet displacement contribution to the stress equation. zero body force.
%
% Neumann traction constraints are handled separately.

    f_sigma = zeros(nStressDofs, 1);
    f_u = zeros(nDispDofs, 1);

    for f = 1:numel(face_struct)

        face = face_struct(f);

        if ~face.is_boundary
            continue;
        end

        mask = face.displacement_mask(:);

        if ~any(mask)
            continue;
        end

        geom = face_global_geom(f);

        xf = face.center(:);
        Af = face.area;

        qp = face.quad_points;
        w = face.quad_weights;

        stressDofs = face.stress_dofs;

        % global normal relative to exterior outward normal
        e = face.cells(face.cells > 0);
        assert(numel(e) == 1, 'Boundary face %d must have one adjacent cell.', f);

        outwardDirection = xf - cell_struct(e).center(:);
        normalSign = sign(dot(geom.n, outwardDirection));

        localRHS = zeros(6,1);

        for q = 1:numel(w)

            xq = qp(:,q);

            uD = face.displacement_value(:);
            uD(~mask) = 0;

            Phi = evaluateFaceBasis(xq, xf, geom.Qf, geom.t1, geom.t2, geom.n, Af, geom.m20, geom.m02);

            localRHS = localRHS + w(q)*normalSign*(Phi.'*uD);

        end

        f_sigma(stressDofs) = f_sigma(stressDofs) + localRHS;

    end

    RHS = [f_sigma;
           f_u;
           zeros(nContactTractionDofs,1);
           zeros(nInterfaceDispDofs,1)];

end

