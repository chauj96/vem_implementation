function K_stab = computeStabilization(cell_struct, face_struct, B_local, P_local)
% construct the local stabilization matrix K_E^{stab}
%
% k_E : infinity norm of the cell compliance C^{-1}
% h_E : diameter of cell E
%
% scaling factor: k_E * h_E

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
    
        % projected polynomial stresses
        proj = projToMatrices(P_local{e});

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
            tractionProj = zeros(3,nDofs);
    
            for j = 1:nDofs
                tractionProj(:,j) = proj{j} * n;
            end
    
            % six DOFs belonging to current face
            cols = 6*(lf-1) + (1:6);
    
    
            for q = 1:numel(w)
    
                xq = qp(:,q);
    
                % start with - projected traction
                residual = -tractionProj;
    
                % all six actual face basis functions simultaneously
                Phi = evaluateFaceBasis( ...
                    xq, xf, ...
                    Qf, t1, t2, n, ...
                    Af, m20, m02);
    
                residual(:,cols) = residual(:,cols) + Phi;
    
                % complete stabilization contribution
                K_E = K_E + w(q) * (residual.' * residual);
    
            end
    
        end
    
        K_stab{e} = kE * hE * K_E;
    
    end

end
