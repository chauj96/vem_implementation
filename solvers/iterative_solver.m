function [sigma_h, u_h] = iterative_solver(A_global, rhs, cell_struct, face_struct)
% solve the mixed-VEM saddle-point system using right-preconditioned GMRES

    % solver parameters
    solverTol = 1e-11;
    maxIter = 1000;
    reorth = 1;
    nStressDofs = 6*numel(face_struct);
    nDispDofs = 6*numel(cell_struct);
    
    % block DOFs
    stressDofs = 1:nStressDofs;
    dispDofs = nStressDofs + (1:nDispDofs);
    
    % extract saddle-point blocks
    A_ss = A_global(stressDofs, stressDofs);
    A_su = A_global(stressDofs, dispDofs);
    A_us = A_global(dispDofs, stressDofs);
    A_uu = A_global(dispDofs, dispDofs);
    
    % approximate Schur complement
    Dinv = spdiags(1 ./ diag(A_ss), 0, nStressDofs,nStressDofs);
    S_approx = A_uu - A_us * Dinv * A_su;
    
    % factorize preconditioner blocks
    F_ss = factorize(A_ss);
    F_S = factorize(S_approx);
    
    % right-preconditioned GMRES
    x0 = zeros(size(rhs));
    
    [sol, flag, total_iters] = gmres_r( ...
        @(v) A_global * v, ...
        rhs, ...
        solverTol, ...
        maxIter, ...
        reorth, ...
        @(v) blockPrecMixedVEM(v,F_ss,A_us,F_S,nStressDofs), ...
        x0);
    
    fprintf('GMRES flag       : %d\n',flag);
    fprintf('GMRES iterations : %d\n',total_iters);
    fprintf('Final residual   : %.3e\n', ...
        norm(A_global*sol-rhs)/norm(rhs));
    
    % extract solution
    sigma_h = sol(1:nStressDofs);
    u_h = sol(nStressDofs+1:end);

end


function y = blockPrecMixedVEM(r,F_ss,A_us,F_S,nStressDofs)

    r_sigma = r(1:nStressDofs);
    r_u = r(nStressDofs+1:end);
    
    y_sigma = F_ss \ r_sigma;
    y_u = F_S \ (r_u - A_us*y_sigma);
    
    y = [y_sigma;
         y_u];

end