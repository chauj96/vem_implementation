function Ac = computeInterfaceAc(fracture, face_struct, face_global_geom, cn, cT)
% assemble interface compliance matrix A_c
%
% a_c(tc, sc) = integral_Gamma sc' * Dc * tc dS
%
% contact traction uses one global/reference face basis per fracture face.
%
% Ac: nContactTractionDofs x nContactTractionDofs

    nTc = fracture.nTractionDofs;

    Ac = sparse(nTc, nTc);

    assert(cn > 0 && cT > 0, 'Contact stiffnesses must be positive.');

    for k = 1:fracture.nFaces

        f = fracture.faces(k);

        % global/reference face geometry
        geom = face_global_geom(f);

        n = geom.n(:);
        n = n / norm(n);

        % normal and tangential compliance
        Pn = n*n.';
        Pt = eye(3) - Pn;

        Dc = (1/cn)*Pn + (1/cT)*Pt;

        xf = face_struct(f).center(:);
        Af = face_struct(f).area;

        qp = face_struct(f).quad_points;
        w = face_struct(f).quad_weights;

        % local 6x6 compliance matrix
        Acf = zeros(6,6);

        for q = 1:numel(w)

            Psi = evaluateFaceBasis(qp(:,q), xf, geom.Qf, geom.t1, geom.t2, geom.n, Af, geom.m20, geom.m02);

            Acf = Acf + w(q)*(Psi.'*Dc*Psi);

        end

        % Global contact traction DOFs
        dofs = fracture.traction_dofs(k,:);

        Ac(dofs,dofs) = Ac(dofs,dofs) + Acf;

    end
end
