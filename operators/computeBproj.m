function Bproj = computeBproj(cell_struct, face_struct, D_local, B_local)
% construct the projection RHS matrix
% Bproj = - volume contribution + face contribution
%
% Both contributions are sums of w_q G_q' X_q over quadrature points. Stacking
% the point blocks vertically turns each into a single product,
%
%   sum_q w_q G_q' X_q = G_stack' * (w_stack .* X_stack)
%
% so the per-point function calls and small matrix products are replaced by one
% batched basis evaluation and one GEMM per cell and per face

    nCells = numel(cell_struct);
    Bproj = cell(nCells,1);

    for e = 1:nCells

        xE = cell_struct(e).center(:);

        qp = cell_struct(e).quad_points;
        w = cell_struct(e).quad_weights;

        D_E = D_local{e};

        nCols = size(D_E,2);

        % VOLUME CONTRIBUTION
        Alpha = D_E(1:3,:);
        Omega = D_E(4:6,:);

        r = qp - xE;

        nq = numel(w);

        % div phi_j = alpha_j + omega_j x r, stacked by quadrature point
        divPhi = zeros(3*nq, nCols);

        divPhi(1:3:end,:) = Alpha(1,:) + r(3,:).'*Omega(2,:) - r(2,:).'*Omega(3,:);
        divPhi(2:3:end,:) = Alpha(2,:) + r(1,:).'*Omega(3,:) - r(3,:).'*Omega(1,:);
        divPhi(3:3:end,:) = Alpha(3,:) + r(2,:).'*Omega(1,:) - r(1,:).'*Omega(2,:);

        G = projectionBasisBatch(qp, xE);

        wStack = repelem(w(:), 3);

        Bproj_E = -(G.' * (wStack .* divPhi));

        % FACE CONTRIBUTION
        nElementFaces = numel(cell_struct(e).faces);

        for lf = 1:nElementFaces

            f = cell_struct(e).faces(lf);

            geom = B_local{e}.geom(lf);

            xf = face_struct(f).center(:);
            Af = face_struct(f).area;

            qpFace = face_struct(f).quad_points;
            wFace = face_struct(f).quad_weights;

            % six columns associated with local face lf
            cols = 6*(lf-1) + (1:6);

            Gf = projectionBasisBatch(qpFace, xE);

            Phi = faceBasisBatch(qpFace, xf, geom.Qf, ...
                                 geom.t1, geom.t2, geom.n, ...
                                 Af, geom.m20, geom.m02);

            wFaceStack = repelem(wFace(:), 3);

            Bproj_E(:,cols) = Bproj_E(:,cols) + Gf.' * (wFaceStack .* Phi);

        end

        Bproj{e}.matrix = Bproj_E;

    end

end
