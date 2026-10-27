function RHS = assembleContactRHS(cell_struct, face_struct, face_global_geom, nStressDofs, nDispDofs, nContactTractionDofs, nInterfaceDispDofs, xFrac)
% assemble RHS for constant-stress contact patch test.
%
% prescribed displacement on the entire exterior boundary.
% zero body force.
%
% unknown ordering:
%   [sigma; u; tc; uc]

    f_sigma = zeros(nStressDofs, 1);
    f_u = zeros(nDispDofs, 1);

    for f = 1:numel(face_struct)

        face = face_struct(f);

        if ~face.is_displacement_bc
            continue;
        end

        % adjacent cell determines minus/plus region
        e = face.cells(1);
        isPlus = cell_struct(e).center(1) > xFrac;

        geom = face_global_geom(f);

        xf = face.center(:);
        Af = face.area;

        qp = face.quad_points;
        w = face.quad_weights;

        % global face stress DOFs
        stressDofs = face.stress_dofs;

        % outward normal relative to global face normal
        outwardDirection = xf - cell_struct(e).center(:);
        normalSign = sign(dot(geom.n, outwardDirection));

        localRHS = zeros(6,1);

        for q = 1:numel(w)

            xq = qp(:,q);

            % exact prescribed displacement
            % uD = [0.4*xq(1);
            %      -0.1*xq(2);
            %      -0.1*xq(3)];
            % 
            % if isPlus
            %     uD(1) = uD(1) + 1;
            % end

            uD = [xq(1)/3 + double(isPlus); 0; 0];

            % global traction basis
            Phi = evaluateFaceBasis(xq, xf, geom.Qf, geom.t1, geom.t2, geom.n, Af, geom.m20, geom.m02);

            % boundary displacement contribution
            localRHS = localRHS + w(q) * normalSign * (Phi.' * uD);

        end

        f_sigma(stressDofs) = f_sigma(stressDofs) + localRHS;

    end

    % no external RHS contributions for contact equations
    RHS = [f_sigma;
           f_u;
           zeros(nContactTractionDofs,1);
           zeros(nInterfaceDispDofs,1)];

end
