function Bproj = computeBproj(cell_struct, face_struct, D_local, B_local)
% Construct the projection RHS matrix
% Bproj = - volume contribution + face contribution

    nCells = numel(cell_struct);
    Bproj = cell(nCells,1);

    for e = 1:nCells

        xE = cell_struct(e).center(:);

        qp = cell_struct(e).quad_points;
        w  = cell_struct(e).quad_weights;

        D_E = D_local{e};

        nCols = size(D_E,2);

        Bproj_E = zeros(6,nCols);

        for j = 1:nCols

            % Volume contribution
            alpha = D_E(1:3,j);
            omega = D_E(4:6,j);

            for q = 1:numel(w)

                xq = qp(:,q);
                divPhi = alpha + cross(omega,xq-xE);

                for i = 1:6

                    g = evaluateProjectionBasis(i,xq,xE);
                    Bproj_E(i,j) = Bproj_E(i,j) - w(q)*dot(divPhi,g);

                end
            end

            % Face contribution
            lf = floor((j-1)/6) + 1;
            jf = mod(j-1,6) + 1;

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

            qpFace = face_struct(f).quad_points;
            wFace = face_struct(f).quad_weights;

            for q = 1:numel(wFace)

                xq = qpFace(:,q);

                phi = evaluateFaceBasis(jf,xq,xf,Qf,t1,t2,n,Af,m20,m02,m11);

                for i = 1:6

                    g = evaluateProjectionBasis(i,xq,xE);
                    Bproj_E(i,j) = Bproj_E(i,j) + wFace(q) * dot(phi,g);

                end

            end

        end

        Bproj{e}.matrix = Bproj_E;

    end

end

function g = evaluateProjectionBasis(i, x, xE)
% evaluate the polynomial vector basis g_i satisfying
% epsilon(g_i) = pi_i

    r = x - xE;

    xr = r(1);
    yr = r(2);
    zr = r(3);

    switch i

        case 1
            g = [xr; 0; 0];

        case 2
            g = [0; yr; 0];

        case 3
            g = [0; 0; zr];

        case 4
            g = (1/sqrt(2))*[yr; xr; 0];

        case 5
            g = (1/sqrt(2))*[zr; 0; xr];

        case 6
            g = (1/sqrt(2))*[0; zr; yr];

    end

end