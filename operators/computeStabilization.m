function K_stab = computeStabilization(cell_struct, face_struct, B_local, P_local)
% construct the local stabilization matrix K_E^{stab}
%
% k_E : 1/(2 mu) or infinity norm of C^{-1}
% h_E : diameter of cell E
%
% scaling factor: k_E * h_E
%
% face residual at quadrature point q:  R_q = -N_f P + E_q,
% E_q = Phi_q on the six columns of f and zero elsewhere. Then
%
%   sum_q w_q R_q' R_q = W_f P' N_f' N_f P - Z_f - Z_f' + Q_f
%
%   W_f = sum_q w_q,   Z_f = P' N_f' Phibar_f  on columns of f,
%   Phibar_f = sum_q w_q Phi_q,   Q_f = sum_q w_q Phi_q' Phi_q
%
% summing over f, the P' (sum_f W_f N_f' N_f) P term is formed once per cell

    nCells = numel(cell_struct);
    K_stab = cell(nCells,1);

    mu = 1; % need to come back for scaling!

    s = 1/sqrt(2);

    for e = 1:nCells

        elementFaces = cell_struct(e).faces;
        nFaces = numel(elementFaces);
        nDofs = 6 * nFaces;

        K_E = zeros(nDofs,nDofs);

        kE = 1 / (2*mu);

        % alternative:
        % Cinv = cell_struct(e).Cinv;
        % kE = norm(Cinv,inf);

        hE = cell_struct(e).diameter;

        % projected polynomial stresses
        P = P_local{e};

        % M = sum_f W_f N_f' N_f
        M = zeros(6,6);

        % FACE CONTRIBUTIONS
        for lf = 1:nFaces

            f = elementFaces(lf);

            geom = B_local{e}.geom(lf);

            xf = face_struct(f).center(:);
            Af = face_struct(f).area;

            n = geom.n(:);

            qp = face_struct(f).quad_points;
            w = face_struct(f).quad_weights(:);

            nq = numel(w);

            % N c = S(c) n, the traction of the stress tensor with
            % coefficients c; N P is the projected traction on face f
            N = [n(1),    0,    0,      0, s*n(3), s*n(2); ...
                    0, n(2),    0, s*n(3),      0, s*n(1); ...
                    0,    0, n(3), s*n(2), s*n(1),      0];

            % (3 nq) x 6, blocks Phi_q down the rows
            Phi = evaluateFaceBasisAll(qp, xf, geom.Qf, geom.t1, geom.t2, n, Af, geom.m20, geom.m02);

            % Phibar_f = sum_q w_q Phi_q,   Q_f = sum_q w_q Phi_q' Phi_q
            Phibar = reshape(sum(reshape(Phi, 3, nq, 6) .* reshape(w, 1, nq, 1), 2), 3, 6);
            Qbar = Phi.' * (repelem(w, 3, 1) .* Phi);

            % six DOFs belonging to current face
            cols = 6*(lf-1) + (1:6);

            % Z_f = P' N_f' Phibar_f
            Z = P.' * (N.' * Phibar);

            K_E(:,cols) = K_E(:,cols) - Z;
            K_E(cols,:) = K_E(cols,:) - Z.';
            K_E(cols,cols) = K_E(cols,cols) + Qbar;

            M = M + sum(w) * (N.' * N);

        end

        K_E = K_E + P.' * M * P;

        K_stab{e} = kE * hE * K_E;

    end

end
