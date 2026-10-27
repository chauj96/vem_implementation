function T = computeInterfaceT(fracture, cell_struct, face_struct, B_local, face_global_geom, nStressDofs)
% assemble interface coupling matrix T
%
% t(sigma, vc) = - integral_Gamma (sigma*n) . vc
%
% assumption: interface displacement basis equals the local face traction basis on each side.
%
% T: nInterfaceDispDofs x nStressDofs

    nUc = fracture.nDispDofs;

    T = sparse(nUc, nStressDofs);

    for k = 1:fracture.nFaces

        f = fracture.faces(k);

        cells = [fracture.minus_cells(k), fracture.plus_cells(k)];

        ucDofs = {fracture.disp_minus_dofs(k,:), fracture.disp_plus_dofs(k,:)};

        for side = 1:2

            e = cells(side);

            lf = find(cell_struct(e).faces == f);
            assert(isscalar(lf), 'Fracture face not found in adjacent cell.');

            % local-to-global stress mapping
            colsLocal = 6*(lf-1) + (1:6);
            sigmaDofs = cell_struct(e).stress_dofs(colsLocal);

            % cell-oriented face geometry
            geom = B_local{e}.geom(lf);

            xf = face_struct(f).center(:);
            Af = face_struct(f).area;

            qp = face_struct(f).quad_points;
            w = face_struct(f).quad_weights;

            % face Gram matrix
            Mf = zeros(6,6);

            for q = 1:numel(w)

                Phi = evaluateFaceBasis(qp(:,q), xf, geom.Qf, geom.t1, geom.t2, geom.n, Af, geom.m20, geom.m02);

                Mf = Mf + w(q)*(Phi.'*Phi);

            end

            rows = ucDofs{side};

            Tf = computeFaceTransformation(f, geom, face_global_geom(f), face_struct);

            T(rows, sigmaDofs) = T(rows, sigmaDofs) - Mf * Tf;

        end
    end
end
