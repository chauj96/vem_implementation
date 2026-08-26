function [sigma_h, u_h] = solveMixedVEM(cell_struct, face_struct, face_global_geom, B_local, K_cons, K_stab, rhs_D, rhs_f, lambda, mu, solverOption)
% solve the mixed VEM system

    % ASSEMBLY
    fprintf('\n[# ASSEMBLE GLOBAL SYSTEM]\n');
    
    rhs = [rhs_D;
            rhs_f];
    
    [K_cons_global, K_stab_global, B_global] = assembleGlobalMatrices(cell_struct, face_struct, face_global_geom, B_local, K_cons, K_stab);
    
    K_global = K_cons_global + K_stab_global;
    A_global = [K_global, B_global';
                B_global, sparse(6*numel(cell_struct), 6*numel(cell_struct))];
    
    [A_global, rhs, ~, ~] = assembleNeumann(face_struct, A_global, rhs, face_global_geom, lambda, mu);
    
    % SOLVE: we have two options (direct and iterative)
    fprintf('\n[# SOLVE THE LINEAR SYSTEM OF MIXED VEM: %d x %d]\n', size(A_global,1), size(A_global,2));
    
    switch solverOption
    
        case 'direct'
            sol = A_global \ rhs;
            sigma_h = sol(1:nStressDofs);
            u_h = sol(nStressDofs+1:end);
    
        case 'iterative'
            [sigma_h, u_h] = iterative_solver(A_global, rhs, cell_struct, face_struct);
    
        otherwise
            error('Unknown solver option: %s', solverOption);
    
    end

end