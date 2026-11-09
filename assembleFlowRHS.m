function b = assembleFlowRHS( ...
    cell_struct, fracture, flowDofs, face_struct, exactFun)

    nW  = flowDofs.nW;
    nP  = flowDofs.nP;
    nWF = flowDofs.nWF;
    nPF = flowDofs.nPF;

    bW  = zeros(nW, 1);
    bP  = zeros(nP, 1);
    bWF = zeros(nWF, 1);
    bPF = zeros(nPF, 1);

    %% BULK SOURCE
    assert(numel(cell_struct) == nP, ...
        'Check bulk pressure DOF mapping.');
    
    % Two-point Gauss-Legendre quadrature
    g = [-1, 1] / sqrt(3);
    
    for e = 1:numel(cell_struct)
    
        xc = cell_struct(e).center(:);
        Ve = cell_struct(e).volume;
    
        % Cartesian cubic mesh
        h = Ve^(1/3);
    
        integralF = 0;
    
        for i = 1:2
            for j = 1:2
                for k = 1:2
    
                    xq = xc + (h/2)*[g(i); g(j); g(k)];
    
                    exact = exactFun(xq);
    
                    integralF = integralF + exact.fBulk;
    
                end
            end
        end
    
        bP(e) = (Ve/8) * integralF;
    end


    %% FRACTURE SOURCE

    for k = 1:fracture.nFaces

        f = fracture.faces(k);
        dof = fracture.pressure_dofs(k);

        xf = face_struct(f).center(:);
        Af = fracture.face_areas(k);

        exact = exactFun(xf);

        bPF(dof) = bPF(dof) + Af * exact.fFrac;
    end

    %% GLOBAL RHS

    b = [bW; bP; bWF; bPF];

    fprintf('\n[# ASSEMBLE FLOW RHS: FULL FRACTURE MMS]\n');
    fprintf('  RHS size: %d x 1\n', numel(b));
    fprintf('  Bulk source norm: %.6e\n', norm(bP));
    fprintf('  Fracture source norm: %.6e\n', norm(bPF));
    fprintf('  Status: ASSEMBLED\n');
end


% function b = assembleFlowRHS(cell_struct, fracture, flowDofs)
% 
%     nW = flowDofs.nW;
%     nP = flowDofs.nP;
%     nWF = flowDofs.nWF;
%     nPF = flowDofs.nPF;
% 
%     %% BULK SOURCE: f_m = 0
% 
%     bW = zeros(nW, 1);
%     bP = zeros(nP, 1);
% 
%     %% FRACTURE SOURCE: f_F = -2
% 
%     bWF = zeros(nWF, 1);
%     bPF = zeros(nPF, 1);
% 
%     for k = 1:fracture.nFaces
% 
%         dof = fracture.pressure_dofs(k);
%         Af = fracture.face_areas(k);
% 
%         bPF(dof) = -2 * Af;
%     end
% 
%     %% GLOBAL RHS
% 
%     b = [bW; bP; bWF; bPF];
% 
%     fprintf('\n[# ASSEMBLE FLOW RHS: NONZERO TRANSFER PATCH TEST]\n');
%     fprintf('  RHS size: %d x 1\n', numel(b));
%     fprintf('  Bulk source norm: %.6e\n', norm(bP));
%     fprintf('  Fracture source norm: %.6e\n', norm(bPF));
%     fprintf('  Status: PASS\n');
% end

% function b = assembleFlowRHS(cell_struct, fracture, flowDofs)
% 
%     nW = flowDofs.nW;
%     nP = flowDofs.nP;
%     nWF = flowDofs.nWF;
%     nPF = flowDofs.nPF;
% 
%     bW = zeros(nW, 1);
%     bP = zeros(nP, 1);
%     bWF = zeros(nWF, 1);
%     bPF = zeros(nPF, 1);
% 
%     n = 1.5;
%     q = n + 1;
% 
%     %% BULK SOURCE
% 
%     for e = 1:numel(cell_struct)
% 
%         x = cell_struct(e).center(:);
% 
%         a = x(1) - 0.5;
%         y = x(2);
%         z = x(3);
% 
%         by = y - min(max(y, 0.25), 0.75);
%         bz = z - min(max(z, 0.25), 0.75);
% 
%         d = sqrt(a^2 + by^2 + bz^2);
% 
%         m = 1 + (by ~= 0) + (bz ~= 0);
% 
%         if d == 0
%             fm = 0;
%         else
%             fm = -q * (q + m - 2) * d^(q - 2);
%         end
% 
%         inside = y >= 0.25 && y <= 0.75 ...
%             && z >= 0.25 && z <= 0.75;
% 
%         if inside && d > 0
% 
%             beta1 = y - 0.25;
%             beta2 = y - 0.75;
%             gamma1 = z - 0.25;
%             gamma2 = z - 0.75;
% 
%             omega = beta1^2 * beta2^2 * gamma1^2 * gamma2^2;
% 
%             omegaYY = 2 * (beta1^2 + 4 * beta1 * beta2 + beta2^2) ...
%                 * gamma1^2 * gamma2^2;
% 
%             omegaZZ = 2 * (gamma1^2 + 4 * gamma1 * gamma2 + gamma2^2) ...
%                 * beta1^2 * beta2^2;
% 
%             fm = fm - d * (omegaYY + omegaZZ);
%         end
% 
%         pDof = cell_struct(e).pressure_dof;
% 
%         bP(pDof) = cell_struct(e).volume * fm;
%     end
% 
%     %% FRACTURE SOURCE (3x3 GAUSS QUADRATURE)
% 
%     xi = [-sqrt(3/5), 0, sqrt(3/5)];
%     wi = [5/9, 8/9, 5/9];
% 
%     for k = 1:fracture.nFaces
% 
%         x = fracture.face_centers(k, :);
%         h = sqrt(fracture.face_areas(k));
% 
%         integralF = 0;
% 
%         for i = 1:3
%             for j = 1:3
% 
%                 y = x(2) + 0.5 * h * xi(i);
%                 z = x(3) + 0.5 * h * xi(j);
% 
%                 beta1 = y - 0.25;
%                 beta2 = y - 0.75;
%                 gamma1 = z - 0.25;
%                 gamma2 = z - 0.75;
% 
%                 omega = beta1^2 * beta2^2 * gamma1^2 * gamma2^2;
% 
%                 omegaYY = 2 * (beta1^2 + 4 * beta1 * beta2 + beta2^2) ...
%                     * gamma1^2 * gamma2^2;
% 
%                 omegaZZ = 2 * (gamma1^2 + 4 * gamma1 * gamma2 + gamma2^2) ...
%                     * beta1^2 * beta2^2;
% 
%                 fF = omegaYY + omegaZZ - 2 * omega;
% 
%                 integralF = integralF + wi(i) * wi(j) * fF;
%             end
%         end
% 
%         pFDof = fracture.pressure_dofs(k);
%         bPF(pFDof) = 0.25 * h^2 * integralF;
%     end
% 
%     %% GLOBAL RHS
% 
%     b = [bW; bP; bWF; bPF];
% 
%     fprintf('\n[# ASSEMBLE FLOW RHS]\n');
%     fprintf('  RHS size: %d x 1\n', numel(b));
%     fprintf('  Bulk source norm: %.6e\n', norm(bP));
%     fprintf('  Fracture source norm: %.6e\n', norm(bPF));
%     fprintf('  Status: PASS\n');
% end
