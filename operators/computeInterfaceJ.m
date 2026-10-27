function J = computeInterfaceJ(fracture, cell_struct, face_struct, B_local, face_global_geom)

    J = sparse(fracture.nTractionDofs, fracture.nDispDofs);

    for k = 1:fracture.nFaces

        f = fracture.faces(k);

        xf = face_struct(f).center(:);
        Af = face_struct(f).area;

        qp = face_struct(f).quad_points;
        w = face_struct(f).quad_weights;

        % common contact traction basis
        gc = face_global_geom(f);

        tcDofs = fracture.traction_dofs(k,:);

        cells = [fracture.minus_cells(k), fracture.plus_cells(k)];

        ucDofs = {fracture.disp_minus_dofs(k,:), fracture.disp_plus_dofs(k,:)};

        % J = [ + crossGram_minus, - crossGram_plus ]
        signs = [1, -1];

        for side = 1:2

            e = cells(side);

            lf = find(cell_struct(e).faces == f);

            geom = B_local{e}.geom(lf);

            Jf = zeros(6,6);

            for q = 1:numel(w)

                Psi = evaluateFaceBasis(qp(:,q), xf, gc.Qf, gc.t1, gc.t2, gc.n, Af, gc.m20, gc.m02);
                Phi = evaluateFaceBasis(qp(:,q), xf, geom.Qf, geom.t1, geom.t2, geom.n, Af, geom.m20, geom.m02);

                Jf = Jf + w(q)*(Psi.'*Phi);

            end

            J(tcDofs,ucDofs{side}) = J(tcDofs,ucDofs{side}) + signs(side)*Jf;

        end
    end
end
