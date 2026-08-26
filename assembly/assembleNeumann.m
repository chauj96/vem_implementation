function [A_global, rhs, neuDofs, neuVals] = assembleNeumann(face_struct, A_global, rhs, face_global_geom, lambda, mu)
% apply prescribed Neumann stress DOFs to the mixed global system

    [neuDofs, neuVals] = computeNeumannData(face_struct, face_global_geom, lambda, mu);
    
    % move known Neumann contributions to RHS
    rhs = rhs - A_global(:,neuDofs) * neuVals;
    
    % zero rows and columns of prescribed stress DOFs
    A_global(:,neuDofs) = 0;
    A_global(neuDofs,:) = 0;
    
    % impose sigma_i = prescribed value
    A_global(neuDofs,neuDofs) = speye(numel(neuDofs));
    rhs(neuDofs) = neuVals;

end