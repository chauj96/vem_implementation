function errors = computeSolutionError(cell_struct, face_struct, sigma_h, u_h, face_global_geom, lambda, mu)
% compute stress and displacement DOF errors for the patch test
%
%   u_exact = [x + 1; 0; 0]
%
% using the global face convention for stress DOFs and
% the cell-wise RM(E) representation for displacement DOFs

    nFaces = numel(face_struct);
    nCells = numel(cell_struct);
    
    nStressDofs = 6 * nFaces;
    nDispDofs   = 6 * nCells;
    
    %% EXACT STRESS  
    % u = [x+1; 0; 0]
    %
    % epsilon(u) = diag(1,0,0)
    %
    % sigma = 2*mu*epsilon + lambda*tr(epsilon)*I
    
    sigma_exact = [2*mu + lambda, 0,      0;
                   0,             lambda, 0;
                   0,             0,      lambda];
    
    sigma_exact_global = zeros(nStressDofs, 1);
    
    for f = 1:nFaces
    
        Af = face_struct(f).area;
        geom = face_global_geom(f);
    
        % exact traction in the global/reference face convention
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
    
    %% EXACT DISPLACEMENT
    u_exact = zeros(nDispDofs, 1);
    
    for e = 1:nCells
    
        xE = cell_struct(e).center(:);
    
        % RM(E) representation of u = [x+1; 0; 0]
        coeff = [ ...
            xE(1) + 1;
            0;
            0;
            0;
            0;
            0];
    
        idx = 6*(e-1) + (1:6);
        u_exact(idx) = coeff;
    
    end
    
    %% STRESS ERROR
    
    stressDiff = sigma_h - sigma_exact_global;
    
    errors.stress.abs = norm(stressDiff);
    errors.stress.rel = norm(stressDiff) / norm(sigma_exact_global);
    
    %% DISPLACEMENT ERROR
    
    dispDiff = u_h - u_exact;
    
    errors.displacement.abs = norm(dispDiff);
    errors.displacement.rel = norm(dispDiff) / norm(u_exact);
    
    %% STORE EXACT SOLUTIONS
    
    errors.sigma_exact = sigma_exact;
    errors.sigma_exact_global = sigma_exact_global;
    errors.u_exact = u_exact;
    
    %% PRINT RESULTS
    
    fprintf('\nSOLUTION ERROR CHECK\n');
    fprintf('=====================================================\n');
    
    fprintf('Stress absolute error      : %.3e\n', errors.stress.abs);
    fprintf('Stress relative error      : %.3e\n', errors.stress.rel);
    
    fprintf('\n');
    
    fprintf('Displacement absolute error: %.3e\n', errors.displacement.abs);
    fprintf('Displacement relative error: %.3e\n', errors.displacement.rel);
    
    fprintf('=====================================================\n');

end