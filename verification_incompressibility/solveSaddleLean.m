function [sol, info] = solveSaddleLean(A, rhs, nS, nD, tol, maxit)
% memory-lean solve of the mixed saddle system
%
%   [ K  B'] [sigma]   [g]
%   [ B  0 ] [  u  ] = [f]
%
% symmetric two-stage diagonal scaling, then MINRES with a block-diagonal
% incomplete-Cholesky preconditioner. MINRES has a three-term recurrence, so
% its memory is O(n). the alternatives are much heavier at the 3D sizes:
% a sparse factorization pays the fill-in of a 1e5-DOF saddle system, and
% solvers/gmres_r preallocates an n-by-maxit Krylov basis (823 MB at 1e5 DOFs)
% on top of an exact factorization of K.
%
% info.relres is measured on the ORIGINAL system, not the scaled one

    if nargin < 5 || isempty(tol),   tol   = 1e-11; end
    if nargin < 6 || isempty(maxit), maxit = 2000;  end

    K = A(1:nS, 1:nS);
    B = A(nS+1:end, 1:nS);

    % stage 1: unit diagonal on the stress block
    d = full(diag(K));
    d(d <= 0) = 1;
    s1 = 1 ./ sqrt(d);

    Ks = (s1 .* K) .* s1.';
    Bs = B .* s1.';

    % stage 2: unit diagonal on the approximate Schur complement
    Shat = Bs * Bs.';
    d2 = full(diag(Shat));
    d2(d2 <= 0) = 1;
    s2 = 1 ./ sqrt(d2);

    Bs = s2 .* Bs;
    Shat = (s2 .* Shat) .* s2.';

    clear K B

    As = [Ks, Bs.'; Bs, sparse(nD, nD)];
    bs = [s1 .* rhs(1:nS); s2 .* rhs(nS+1:end)];

    LK = robustIchol(Ks);
    LS = robustIchol(Shat);

    clear Shat

    prec = @(r) [LK.' \ (LK \ r(1:nS)); LS.' \ (LS \ r(nS+1:end))];

    [y, flag, relres, iter] = minres(As, bs, tol, maxit, prec);

    clear As

    sol = [s1 .* y(1:nS); s2 .* y(nS+1:end)];

    info.flag = flag;
    info.iterations = iter;
    info.scaledRelres = relres;
    info.relres = norm(rhs - A*sol) / norm(rhs);

end

function L = robustIchol(M)
% incomplete Cholesky, retried with a growing diagonal shift on breakdown

    shift = 0;

    for attempt = 1:8
        try
            L = ichol(M, struct('type','ict','droptol',1e-3,'diagcomp',shift));
            return;
        catch
            if shift == 0
                shift = 1e-4;
            else
                shift = shift * 10;
            end
        end
    end

    % last resort: Jacobi
    L = spdiags(sqrt(max(full(diag(M)), eps)), 0, size(M,1), size(M,2));

end
