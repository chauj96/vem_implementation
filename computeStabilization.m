function K_stab = computeStabilization(cell_struct, face_struct, B_local, P_local)
% construct the local stabilization matrix K_E^{stab}
% k_E : 1/(2 mu) or infinity norm of C^{-1}
% h_E : diameter of cell E

    nCells = numel(cell_struct);
    K_stab = cell(nCells, 1);

    mu = 1; % need to come back for scaling!

    for e = 1:nCells
        nFaces = numel(cell_struct(e).faces);
        nDofs = 6 * nFaces;

        K_E = zeros(nDofs, nDofs);
        
        kE = 1 / (2*mu);
        hE = cell_struct(e).diameter;

        % projected polynomial stresses
        proj = projToMatrices(P_local{e});

        for lf = 1:nFaces

            f = cell_struct(e).faces(lf);
            geom = B_local{e}.geom(lf);

            xf = face_struct(f).center(:);
            Af = face_struct(f).area;

            t1 = geom.t1;
            t2 = geom.t2;
            n = geom.n;
            Qf = geom.Qf;

            m20 = geom.m20;
            m02 = geom.m02;
            m11 = geom.m11;

            qp = face_struct(f).quad_points;
            w = face_struct(f).quad_weights;

            for q = 1:numel(w)

                xq = qp(:,q);
                residual = zeros(3, nDofs);

                for j = 1:nDofs

                    face_j = floor((j-1)/6) + 1;
                    local_j = mod(j-1,6) + 1;

                    if face_j == lf
                        phi = evaluateFaceBasis(local_j, xq, xf, Qf, t1, t2, n, Af, m20, m02, m11);
                    else
                        phi = zeros(3,1);
                    end

                    % projected traction
                    tractionProj = proj{j} * n;
                    residual(:,j) = phi - tractionProj;

                end

                K_E = K_E + w(q) * (residual.' * residual);
            end
        end

        K_E = kE * hE * K_E;
        K_stab{e} = K_E;

    end
end