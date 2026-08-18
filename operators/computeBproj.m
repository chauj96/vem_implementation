function Bproj = computeBproj(cell_struct, face_struct, D_local, B_local)
% construct the projection RHS matrix
% Bproj = - volume contribution + face contribution

    nCells = numel(cell_struct);
    Bproj = cell(nCells,1);
    
    for e = 1:nCells
    
        xE = cell_struct(e).center(:);
    
        qp = cell_struct(e).quad_points;
        w = cell_struct(e).quad_weights;
    
        D_E = D_local{e};
    
        nCols = size(D_E,2);
    
        Bproj_E = zeros(6,nCols);
    
        % VOLUME CONTRIBUTION
        Alpha = D_E(1:3,:);
        Omega = D_E(4:6,:);
    
        for q = 1:numel(w)
    
            xq = qp(:,q);
            r = xq - xE;
    
            G = evaluateProjectionBasis(xq, xE);
    
            % omega_j x r for every column j
            crossOmegaR = [ ...
                Omega(2,:)*r(3) - Omega(3,:)*r(2);
                Omega(3,:)*r(1) - Omega(1,:)*r(3);
                Omega(1,:)*r(2) - Omega(2,:)*r(1)];
    
            divPhi = Alpha + crossOmegaR;
    
            % all columns j simultaneously
            Bproj_E = Bproj_E - w(q) * (G' * divPhi);
    
        end
    
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
    
            % complete 6-by-6 contribution from this face
            faceBlock = zeros(6,6);
    
            for q = 1:numel(wFace)
    
                xq = qpFace(:,q);
    
                % all six projection basis functions: G = [g1 ... g6] -> 3-by-6
                G = evaluateProjectionBasis(xq, xE);
    
                % all six face traction basis functions: Phi = [phi1 ... phi6] -> 3-by-6
                Phi = evaluateFaceBasis( ...
                    xq, xf, ...
                    geom.Qf, ...
                    geom.t1, geom.t2, geom.n, ...
                    Af, ...
                    geom.m20, geom.m02);
    
                % G' * Phi is the complete 6-by-6 face block
                faceBlock = faceBlock + wFace(q) * (G' * Phi);
    
            end
    
            Bproj_E(:,cols) = Bproj_E(:,cols) + faceBlock;
    
        end
    
        Bproj{e}.matrix = Bproj_E;
    
    end

end

function G = evaluateProjectionBasis(x, xE)
% evaluate all six polynomial vector basis functions simultaneously.

    r = x - xE;
    
    xr = r(1);
    yr = r(2);
    zr = r(3);
    
    s = 1/sqrt(2);
    
    % shear ordering must match projToMatrices: 4 -> yz, 5 -> xz, 6 -> xy.
    % columns 4 and 6 were previously swapped, so shear stress coefficients were
    % produced as xy/yz but interpreted as yz/xy. invisible for diagonal stress,
    % an O(1) error for any shear
    G = [ ...
        xr,  0,   0,   0,    s*zr, s*yr;
        0,   yr,  0,   s*zr, 0,    s*xr;
        0,   0,   zr,  s*yr, s*xr, 0    ];

end