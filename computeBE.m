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


% helper functions 1:
function [t1, t2, Qf] = constructFaceFrame(X, n)
% construct an orthonormal local frame {t1,t2,n}.

    nVertices = size(X, 1);
    t1 = []; 

    for k = 1:nVertices

        kp1 = mod(k, nVertices) + 1;
        edge = X(kp1,:)' - X(k,:)';

        edge = edge - dot(edge, n) * n;

        if norm(edge) > 1e-14
            t1 = edge / norm(edge);
            break;
        end
    end

    if isempty(t1)
        error('Could not construct t1: all face edges are degenerate.');
    end

    t2 = cross(n, t1);
    t2 = t2 / norm(t2);

    % reorthogonalize t1 to reduce floating-point drift
    t1 = cross(t2, n);
    t1 = t1 / norm(t1);

    Qf = [t1, t2, n];
end

% helper function 2:
function [m20, m02, m11] = computeFaceMoments(X, xf, Qf)
% compute the second-order face moments

    nVertices = size(X, 1);

    % local coordinates of all face vertices
    Xshift = X' - xf;
    Xlocal = Qf' * Xshift;

    xtilde = Xlocal(1,:);
    ytilde = Xlocal(2,:);

    % planarity check
    normalCoordinates = Xlocal(3,:);

    faceScale = max(vecnorm(Xshift, 2, 1));

    if max(abs(normalCoordinates)) > 1e-10 * max(faceScale, eps)
        warning('Face vertices are not exactly coplanar.');
    end

    m20 = 0.0;
    m02 = 0.0;
    m11 = 0.0;

    for k = 1:nVertices

        kp1 = mod(k, nVertices) + 1;

        x1 = xtilde(k);
        y1 = ytilde(k);

        x2 = xtilde(kp1);
        y2 = ytilde(kp1);

        Jk = abs(x1*y2 - x2*y1);

        m20 = m20 + Jk * (x1^2 + x1*x2 + x2^2)/12.0;
        m02 = m02 + Jk * (y1^2 + y1*y2 + y2^2)/12.0;
        m11 = m11 + Jk * (x1*y1/12.0 + x1*y2/24.0 + x2*y1/24.0 + x2*y2/12.0);
    end
end