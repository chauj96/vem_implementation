function B_local = computeBE(cell_struct, face_struct, V3)

    % construct the local divergence bilinear matrix B_E for every cell
    
    nCells = length(cell_struct);
    B_local = cell(nCells, 1);
    
    for e = 1:nCells
        xE = cell_struct(e).center(:); % barycenter of cell E
        
        elementFaces = cell_struct(e).faces;
        nElementFaces = numel(elementFaces);
    
        B_E = zeros(6, 6 * nElementFaces);
    
        for lf = 1:nElementFaces
            f = elementFaces(lf);
    
            xf = face_struct(f).center(:);
            Af = face_struct(f).area;
    
            % current element's outward unit normal
            n = cell_struct(e).face_normals(lf, :)';
            n = n / norm(n);
    
            vertexIds = face_struct(f).verts(:);
            X = V3(vertexIds, :);  % number of face vertices x 3
    
            % construct local orthonormal face frame
            [t1, t2, Qf] = constructFaceFrame(X, n);
    
            % compute second order face moments
            [m20, m02, m11] = computeFaceMoments(X, xf, Qf);
    
            % compute rotational rhs vectors
            df = xf - xE;
            
            b1 = cross(df, t1);
            b2 = cross(df, t2);
            b3 = cross(df, n);
    
            b4 = ((m20 + m02) / Af^3) * n;
            b5 = (m11*t1 - m20*t2) / Af^3;
            b6 = (m02*t1 - m11*t2) / Af^3;
    
            Bomega_f = [b1, b2, b3, b4, b5, b6];
    
            % assemble the 6 by 6 face block
            Bface = [t1, t2, n, zeros(3,3); Bomega_f];
    
            cols = 6 * (lf - 1) + (1:6);
            B_E(:, cols) = Bface;
        end
    
        B_local{e} = B_E;
    end
end
