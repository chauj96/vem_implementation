function problem = problemPatchTest(lambda, mu)
% linear patch test: u = [x+1; 0; 0], constant stress, zero body force
%
% used to validate the convergence driver itself: every level must reproduce
% this solution to machine precision, independent of the mesh

    if nargin < 1 || isempty(lambda), lambda = 1; end
    if nargin < 2 || isempty(mu), mu = 1; end

    S0 = diag([2*mu + lambda, lambda, lambda]);

    problem.dim = 3;
    problem.name = sprintf('patch test (linear u), lambda=%g mu=%g', lambda, mu);
    problem.mu     = @(x) mu * ones(1, size(x,2));
    problem.lambda = @(x) lambda * ones(1, size(x,2));
    problem.u      = @(x) [x(1,:) + 1; zeros(2, size(x,2))];
    problem.sigma  = @(x) repmat(S0, 1, 1, size(x,2));
    problem.f      = @(x) zeros(3, size(x,2));

end
