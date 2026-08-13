function [A_global, rhs, neuDofs, neuVals] = assembleNeumann(face_struct, A_global, rhs, face_global_geom, lambda, mu)
% apply prescribed Neumann stress DOFs to the global system
% for the exact displacement u = [x+1; 0; 0]
    
    % exact stress corresponding to u = [x+1; 0; 0]
    sigma = [2*mu + lambda, 0,      0;
             0,             lambda, 0;
             0,             0,      lambda];
    
    neuDofs = [];
    neuVals = [];
    
    for f = 1:numel(face_struct)
    
        if ~face_struct(f).is_neumann
            continue;
        end
    
        Af = face_struct(f).area;
    
        % global/reference face geometry
        geom = face_global_geom(f);
    
        ng = geom.n;
        t1g = geom.t1;
        t2g = geom.t2;
    
        % traction in the global face convention
        traction = sigma * ng;
    
        % coefficients in the global face traction basis
        coeff = [ ...
            Af * dot(traction, t1g);
            Af * dot(traction, t2g);
            Af * dot(traction, ng);
            0;
            0;
            0];
    
        % corresponding global stress DOFs
        idx = 6*(f-1) + (1:6);
    
        neuDofs = [neuDofs, idx];
        neuVals = [neuVals; coeff];
    
    end
    
    % move known Neumann contributions to RHS
    rhs = rhs - A_global(:, neuDofs) * neuVals;
    
    % zero rows and columns of prescribed stress DOFs
    A_global(:, neuDofs) = 0;
    A_global(neuDofs, :) = 0;
    
    % impose sigma_i = prescribed value
    A_global(neuDofs, neuDofs) = speye(numel(neuDofs));
    rhs(neuDofs) = neuVals;

end