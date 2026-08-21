function K_stab = computeStabilization(cell_struct, face_struct, B_local, P_local)
% construct the local stabilization matrix K_E^{stab}
%
% k_E : infinity norm of the cell compliance C^{-1}
% h_E : diameter of cell E
%
% scaling factor: k_E * h_E
%
% On each face the integrand is R_q = -T + S_q, with T the projected traction
% (constant in q) and S_q the face basis Phi_q placed in this face's columns.
% Expanding
%
%   sum_q w_q R_q' R_q = |f| T'T - T'S~ - S~'T + sum_q w_q S_q'S_q
%
% with S~ = sum_q w_q S_q leaves a single nDofs-by-nDofs product per face
% instead of one per quadrature point; the remaining sums are 3-by-6 and 6-by-6

    nCells = numel(cell_struct);
    K_stab = cell(nCells,1);

    for e = 1:nCells

        nFaces = numel(cell_struct(e).faces);
        nDofs = 6 * nFaces;

        K_E = zeros(nDofs,nDofs);

        % k_E = 1/(2 mu), read off the cell's own compliance: Cinv(4,4) = 1/(2 mu)
        % exactly. this is the DEVIATORIC compliance and carries no lambda, so the
        % stabilization stays fixed as the material approaches incompressibility.
        % the normal-block entries a, b both contain lambda, so ||C^-1||_inf or
        % tr(C^-1) would make the stabilization drift with lambda
        kE = cell_struct(e).Cinv(4,4);

        hE = cell_struct(e).diameter;

        % projected polynomial stresses, stacked so that the projected traction
        % for every DOF is one product rather than nDofs small ones
        proj = projToMatrices(P_local{e});
        projStack = vertcat(proj{:});                  % (3 nDofs)-by-3

        % FACE CONTRIBUTIONS
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

            qp = face_struct(f).quad_points;
            w = face_struct(f).quad_weights;

            % projected traction for ALL element stress DOFs
            tractionProj = reshape(projStack * n, 3, nDofs);

            % six DOFs belonging to current face
            cols = 6*(lf-1) + (1:6);

            % face basis moments: M1 = sum_q w_q Phi_q, M2 = sum_q w_q Phi_q'Phi_q
            M1 = zeros(3,6);
            M2 = zeros(6,6);
            wSum = 0;

            for q = 1:numel(w)

                Phi = evaluateFaceBasis( ...
                    qp(:,q), xf, ...
                    Qf, t1, t2, n, ...
                    Af, m20, m02);

                M1 = M1 + w(q) * Phi;
                M2 = M2 + w(q) * (Phi.' * Phi);
                wSum = wSum + w(q);

            end

            TM = tractionProj.' * M1;                  % nDofs-by-6

            K_E = K_E + wSum * (tractionProj.' * tractionProj);

            K_E(:,cols) = K_E(:,cols) - TM;
            K_E(cols,:) = K_E(cols,:) - TM.';
            K_E(cols,cols) = K_E(cols,cols) + M2;

        end

        K_stab{e} = kE * hE * K_E;

    end

end
