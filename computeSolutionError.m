function errors = computeSolutionError(cell_struct, face_struct, sigma_h, u_h, face_global_geom, lambda, mu)
% compute stress and displacement errors for the patch test
%
% exact displacement:
%       u(x) = [x + 1; 0; 0]
%
% exact displacement is compared with its cell-wise L2 projection onto the rigid-motion space RM(E)

    nFaces = numel(face_struct);
    nCells = numel(cell_struct);
    
    nStressDofs = 6 * nFaces;
    nDispDofs = 6 * nCells;
    
    %% exact stress
    sigma_exact = [2*mu + lambda, 0,      0;
                   0,             lambda, 0;
                   0,             0,      lambda];
    
    sigma_exact_global = zeros(nStressDofs, 1);
    
    for f = 1:nFaces
    
        Af = face_struct(f).area;
        geom = face_global_geom(f);
    
        % exact traction in global/reference face convention
        traction = sigma_exact * geom.n;
    
        % exact global stress DOFs
        coeff = [ ...
            Af * dot(traction, geom.t1);
            Af * dot(traction, geom.t2);
            Af * dot(traction, geom.n);
            0;
            0;
            0];
    
        idx = 6*(f-1) + (1:6);
        sigma_exact_global(idx) = coeff;
    
    end
    
    %% exact displacement
    %  L2 projection of u(x) = [x+1; 0; 0] onto RM(E)
    %
    %  RM(E):
    %
    %  R1 = [1;0;0]
    %  R2 = [0;1;0]
    %  R3 = [0;0;1]
    %
    %  R4 = e1 x (x-xE)
    %  R5 = e2 x (x-xE)
    %  R6 = e3 x (x-xE)
    
    u_exact = zeros(nDispDofs, 1);
    
    e1 = [1;0;0];
    e2 = [0;1;0];
    e3 = [0;0;1];
    
    for e = 1:nCells
    
        xE = cell_struct(e).center(:);
    
        % translation coefficients  
        % since xE is the cell barycenter: average_E u(x) = [xE(1)+1; 0; 0]
    
        alpha = [xE(1) + 1;
                 0;
                 0];
    
        % rotation coefficients
        % rotational Gram matrix:
        %
        % Grot(i,j) = integral_E Ri * Rj dV
        %
        % RHS:
        %
        % brot(i) = integral_E u_exact * Ri dV
    
        Grot = zeros(3,3);
        brot = zeros(3,1);
    
        qp = cell_struct(e).quad_points;
        qw = cell_struct(e).quad_weights;
    
        for q = 1:size(qp,2)
    
            xq = qp(:,q);
            wq = qw(q);
    
            r = xq - xE;
    
            % exact displacement at quadrature point
            uq = [xq(1) + 1;
                  0;
                  0];
    
            % rotational RM basis
            R4 = cross(e1, r);
            R5 = cross(e2, r);
            R6 = cross(e3, r);
    
            R = [R4, R5, R6];
    
            % rotational Gram matrix
            Grot = Grot + wq * (R' * R);
    
            % projection RHS
            brot = brot + wq * (R' * uq);
    
        end
    
        % rotational coefficients
        omega = Grot \ brot;
    
        % complete six RM(E) coefficients
        coeff = [alpha;
                 omega];
    
        idx = 6*(e-1) + (1:6);
        u_exact(idx) = coeff;
    
    end
    
    % stress error    
    stressDiff = sigma_h - sigma_exact_global;
    
    errors.stress.abs = norm(stressDiff);
    
    sigmaNorm = norm(sigma_exact_global);
    
    if sigmaNorm > 0
        errors.stress.rel = errors.stress.abs / sigmaNorm;
    else
        errors.stress.rel = NaN;
    end

    %  displacement error  
    dispDiff = u_h - u_exact;
    
    errors.displacement.abs = norm(dispDiff);
    
    uNorm = norm(u_exact);
    
    if uNorm > 0
        errors.displacement.rel = errors.displacement.abs / uNorm;
    else
        errors.displacement.rel = NaN;
    end
    
    errors.sigma_exact = sigma_exact;
    errors.sigma_exact_global = sigma_exact_global;
    errors.u_exact = u_exact;
    
    fprintf('\n')
    fprintf('=====================================================\n');
    fprintf('# SOLUTION ERROR CHECK\n');
    fprintf('=====================================================\n');
    
    fprintf('\nStress absolute error      : %.3e\n', errors.stress.abs);
    fprintf('Stress relative error      : %.3e\n', errors.stress.rel);

    fprintf('\nDisplacement absolute error: %.3e\n', errors.displacement.abs);
    fprintf('Displacement relative error: %.3e\n', errors.displacement.rel);

    fprintf('=====================================================\n');

end