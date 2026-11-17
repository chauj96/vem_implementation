function [K_bc, RHS_bc, bc] = applyContactBoundaryConstraints(K_global, RHS, cell_struct, face_struct, face_global_geom, nStressDofs)
% impose traction boundary conditions and remove in-plane rigid body modes

    nTotal = size(K_global, 1);

    C = sparse(0, nTotal);
    d = zeros(0, 1);

    %% 1. traction boundary conditions

    for f = 1:numel(face_struct)

        face = face_struct(f);

        if ~face.is_boundary
            continue;
        end

        mask = logical(face.traction_mask(:));

        if ~any(mask)
            continue;
        end

        geom = face_global_geom(f);

        % determine outward normal sign
        e = face.cells(face.cells > 0);
        outwardDirection = face.center(:) - cell_struct(e).center(:);
        normalSign = sign(dot(geom.n, outwardDirection));

        assert(normalSign ~= 0, 'Invalid boundary normal.');

        % construct constraints using actual face basis
        qp = face.quad_points;
        w = face.quad_weights;

        M = zeros(6, 6);
        g = zeros(6, 1);

        P = diag(double(mask));
        traction = face.traction_value(:);

        for q = 1:numel(w)

            Phi = evaluateFaceBasis(qp(:,q), face.center(:), geom.Qf, geom.t1, geom.t2, geom.n, face.area, geom.m20, geom.m02);

            M = M + w(q)*(Phi.'*P*Phi);
            g = g + w(q)*(Phi.'*P*traction);

        end

        Cface = normalSign*M;

        % keep independent constraint equations
        [~, pivots] = rref(full(Cface));

        Cnew = sparse(numel(pivots), nTotal);
        Cnew(:, face.stress_dofs(:)) = Cface(pivots,:);

        C = [C; Cnew];
        d = [d; g(pivots)];

    end

    nTractionConstraints = size(C,1);

    %% 2. remove in-plane rigid body modes

    nCells = numel(cell_struct);

    C_RM = sparse(3, nTotal);

    for e = 1:nCells

        xE = cell_struct(e).center(:);

        qp = cell_struct(e).quad_points;
        w = cell_struct(e).quad_weights;

        C_E = zeros(3,6);

        for q = 1:numel(w)

            xq = qp(:,q);

            % evaluate six displacement basis functions
            R = zeros(3,6);

            for i = 1:6
                R(:,i) = evaluateRMBasis(i, xq, xE);
            end

            % global in-plane rigid body modes
            R_global = [1, 0, -xq(2); 0, 1, xq(1); 0, 0, 0];

            % integral constraints
            C_E = C_E + w(q)*(R_global.'*R);

        end

        % global bulk displacement dofs
        dispDofs = nStressDofs + 6*(e-1) + (1:6);

        C_RM(:,dispDofs) = C_E;

    end

    % normalize constraint rows for numerical scaling
    rowNorms = full(sqrt(sum(C_RM.^2,2)));

    assert(all(rowNorms > 0), 'Invalid rigid body constraints.');

    C_RM = spdiags(1./rowNorms,0,3,3)*C_RM;

    % append rigid body constraints
    C = [C; C_RM];
    d = [d; zeros(3,1)];

    %% 3. lagrange multiplier augmentation

    nConstraints = size(C,1);

    K_bc = [K_global, C.'; C, sparse(nConstraints,nConstraints)];
    RHS_bc = [RHS; d];

    bc.nOriginal = nTotal;
    bc.nConstraints = nConstraints;
    bc.nTractionConstraints = nTractionConstraints;
    bc.nRigidBodyConstraints = 3;
    bc.C_RM = C_RM;

    fprintf('\n[# BOUNDARY CONSTRAINTS]\n');
    fprintf('Traction constraints: %d\n', nTractionConstraints);
    fprintf('Rigid body constraints: %d\n', bc.nRigidBodyConstraints);
    fprintf('Total constraints: %d\n', nConstraints);

end

% function [K_bc, RHS_bc, bc] = applyContactBoundaryConstraints(K_global, RHS, cell_struct, face_struct, face_global_geom)
% % impose component-wise constant traction boundary conditions
% 
%     nTotal = size(K_global, 1);
% 
%     C = sparse(0, nTotal);
%     d = zeros(0, 1);
% 
%     for f = 1:numel(face_struct)
% 
%         face = face_struct(f);
% 
%         if ~face.is_boundary
%             continue;
%         end
% 
%         mask = logical(face.traction_mask(:));
% 
%         if ~any(mask)
%             continue;
%         end
% 
%         geom = face_global_geom(f);
% 
%         % determine outward normal sign
%         e = face.cells(face.cells > 0);
%         outwardDirection = face.center(:) - cell_struct(e).center(:);
%         normalSign = sign(dot(geom.n, outwardDirection));
% 
%         assert(normalSign ~= 0, 'Invalid boundary normal.');
% 
%         % construct constraints using the actual face basis
%         qp = face.quad_points;
%         w = face.quad_weights;
% 
%         M = zeros(6, 6);
%         g = zeros(6, 1);
% 
%         P = diag(double(mask));
%         traction = face.traction_value(:);
% 
%         for q = 1:numel(w)
%             Phi = evaluateFaceBasis(qp(:,q), face.center(:), geom.Qf, geom.t1, geom.t2, geom.n, face.area, geom.m20, geom.m02);
%             M = M + w(q) * (Phi.' * P * Phi);
%             g = g + w(q) * (Phi.' * P * traction);
%         end
% 
%         Cface = normalSign * M;
% 
%         % keep independent constraint equations
%         [~, pivots] = rref(full(Cface));
% 
%         Cnew = sparse(numel(pivots), nTotal);
%         Cnew(:, face.stress_dofs(:)) = Cface(pivots, :);
% 
%         C = [C; Cnew];
%         d = [d; g(pivots)];
% 
%     end
% 
%     % enforce constraints using lagrange multipliers
%     nConstraints = size(C, 1);
% 
%     K_bc = [K_global, C.'; C, sparse(nConstraints, nConstraints)];
%     RHS_bc = [RHS; d];
% 
%     bc.nOriginal = nTotal;
%     bc.nConstraints = nConstraints;
% 
% end