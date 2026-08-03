function B_quad = computeBEQuadrature(cell_struct, face_struct, V3)
% construct B_E by quadrature rule

    nCells = numel(cell_struct);
    B_quad = cell(nCells,1);

    for e = 1:nCells

        Bgeom = struct([]);

        xE = cell_struct(e).center(:);

        elementFaces = cell_struct(e).faces;
        nFaces = numel(elementFaces);

        B_E = zeros(6,6*nFaces);

        for lf = 1:nFaces

            f = elementFaces(lf);

            xf = face_struct(f).center(:);
            Af = face_struct(f).area;

            % outward normal
            n = cell_struct(e).face_normals(lf,:)';
            n = n/norm(n);

            % face vertices
            vertexIds = face_struct(f).verts(:);
            X = V3(vertexIds,:);

            % local frame
            [t1,t2,Qf] = constructFaceFrame(X,n);

            % second moments
            [m20,m02,m11] = computeFaceMoments(X,xf,Qf);

            qp = face_struct(f).quad_points;
            w = face_struct(f).quad_weights;

            Bface = zeros(6,6);

            for q = 1:numel(w)

                xq = qp(:,q);

                for j = 1:6

                    phi = evaluateFaceBasis(j,xq,xf,Qf,t1,t2,n,Af,m20,m02,m11);

                    for i = 1:6

                        R = evaluateRMBasis(i,xq,xE);
                        Bface(i,j) = Bface(i,j) + w(q)*dot(phi,R);

                    end
                end
            end

            cols = 6 * (lf-1) + (1:6);
            B_E(:,cols) = Bface;

            % store the same geometric information as computeBE
            Bgeom(lf).face = f;
            Bgeom(lf).t1 = t1;
            Bgeom(lf).t2 = t2;
            Bgeom(lf).n = n;
            Bgeom(lf).Qf = Qf;
            Bgeom(lf).m20 = m20;
            Bgeom(lf).m02 = m02;
            Bgeom(lf).m11 = m11;

        end

        B_quad{e}.matrix = B_E;
        B_quad{e}.geom = Bgeom;

    end

end