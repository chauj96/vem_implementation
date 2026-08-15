function Bproj = computeBproj(cell_struct, face_struct, D_local, B_local)
% construct the projection RHS matrix
% Bproj = - volume contribution + face contribution
%
% with r_q = x_q - x_E and omega x r = -skew(r) omega,
%
%   sum_q w_q G_q' (alpha_j + omega_j x r_q) = Ghat Alpha - Hhat Omega
%
%   Ghat = sum_q w_q G_q',   Hhat = sum_q w_q G_q' skew(r_q),   both 6 x 3
%
% face terms:  sum_q w_q G_q' Phi_q, one product over the stacked points

    nCells = numel(cell_struct);
    Bproj = cell(nCells,1);

    for e = 1:nCells

        xE = cell_struct(e).center(:);

        qp = cell_struct(e).quad_points;
        w = cell_struct(e).quad_weights(:);

        D_E = D_local{e};

        nq = numel(w);

        % VOLUME CONTRIBUTION
        Alpha = D_E(1:3,:);
        Omega = D_E(4:6,:);

        r = qp - xE;

        G = evaluateProjectionBasisAll(qp, xE);

        wG = repelem(w, 3, 1) .* G;

        % skew(r_q) v = r_q x v, blocks stacked down the rows
        S = zeros(3*nq, 3);

        % component rows of each quadrature point block
        i1 = 1:3:3*nq;
        i2 = 2:3:3*nq;
        i3 = 3:3:3*nq;

        S(i1,2) = -r(3,:);  S(i1,3) =  r(2,:);
        S(i2,1) =  r(3,:);  S(i2,3) = -r(1,:);
        S(i3,1) = -r(2,:);  S(i3,2) =  r(1,:);

        Ghat = [sum(wG(i1,:), 1); sum(wG(i2,:), 1); sum(wG(i3,:), 1)].';
        Hhat = wG.' * S;

        Bproj_E = Hhat * Omega - Ghat * Alpha;

        % FACE CONTRIBUTION
        elementFaces = cell_struct(e).faces;
        nElementFaces = numel(elementFaces);

        for lf = 1:nElementFaces

            f = elementFaces(lf);

            geom = B_local{e}.geom(lf);

            xf = face_struct(f).center(:);
            Af = face_struct(f).area;

            qpFace = face_struct(f).quad_points;
            wFace = face_struct(f).quad_weights(:);

            % six columns associated with local face lf
            cols = 6*(lf-1) + (1:6);

            GFace = evaluateProjectionBasisAll(qpFace, xE);

            Phi = evaluateFaceBasisAll( ...
                qpFace, xf, ...
                geom.Qf, ...
                geom.t1, geom.t2, geom.n, ...
                Af, ...
                geom.m20, geom.m02);

            % sum_q w_q G_q' Phi_q
            faceBlock = GFace.' * (repelem(wFace, 3, 1) .* Phi);

            Bproj_E(:,cols) = Bproj_E(:,cols) + faceBlock;

        end

        Bproj{e}.matrix = Bproj_E;

    end

end

function G = evaluateProjectionBasisAll(x, xE)
% evaluate all six polynomial vector basis functions at all quadrature points
%
%   x : 3 x nq quadrature points
%   G : (3 nq) x 6, rows 3(q-1)+(1:3) hold G_q

    nq = size(x, 2);

    r = x - xE(:);

    xr = r(1,:).';
    yr = r(2,:).';
    zr = r(3,:).';

    s = 1/sqrt(2);

    % component rows of each quadrature point block
    i1 = 1:3:3*nq;
    i2 = 2:3:3*nq;
    i3 = 3:3:3*nq;

    G = zeros(3*nq, 6);

    G(i1,1) = xr;
    G(i2,2) = yr;
    G(i3,3) = zr;

    G(i1,4) = s*yr;   G(i2,4) = s*xr;
    G(i1,5) = s*zr;   G(i3,5) = s*xr;
    G(i2,6) = s*zr;   G(i3,6) = s*yr;

end
