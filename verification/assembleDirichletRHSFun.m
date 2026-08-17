function rhs_D = assembleDirichletRHSFun(cell_struct, face_struct, face_global_geom, problem)
% Dirichlet contribution to the stress RHS for a prescribed displacement field
%
%   (rhs_D)_{f,j} = integral_f u_D . (tau_j n_out) dS
%
% variant of assembleDirichletRHS taking u_D from problem instead of a
% hard-coded expression

    nFaces = numel(face_struct);
    rhs_D = zeros(6*nFaces, 1);

    % adjacent cell of each Dirichlet face
    faceCell = zeros(nFaces,1);
    for e = 1:numel(cell_struct)
        for f = cell_struct(e).faces(:).'
            if face_struct(f).is_dirichlet
                faceCell(f) = e;
            end
        end
    end

    dirichletFaces = find([face_struct.is_dirichlet]);

    for k = 1:numel(dirichletFaces)

        f = dirichletFaces(k);

        xf = face_struct(f).center(:);
        Af = face_struct(f).area;

        qp = face_struct(f).quad_points;
        w = face_struct(f).quad_weights(:);

        geom = face_global_geom(f);

        e = faceCell(f);
        assert(e ~= 0, 'Could not find adjacent cell for Dirichlet face %d.', f);

        % orientation of the reference normal relative to the outward normal
        normalSign = sign(dot(geom.n, xf - cell_struct(e).center(:)));
        assert(normalSign ~= 0, 'Could not determine outward orientation for face %d.', f);

        uD = problem.u(qp);                 % 3 x nq

        faceRHS = zeros(6,1);

        for q = 1:numel(w)
            PhiG = evaluateFaceBasis(qp(:,q), xf, geom.Qf, geom.t1, geom.t2, geom.n, ...
                                     Af, geom.m20, geom.m02);
            faceRHS = faceRHS + w(q) * ((normalSign * PhiG).' * uD(:,q));
        end

        rhs_D(6*(f-1) + (1:6)) = faceRHS;

    end

end
