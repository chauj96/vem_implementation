function problem = problemRotatedPatch(lambda, mu)
% the linear patch test with a rigid rotation added
%
%   u = [x + 1 - y; x; 0] = [x+1;0;0]  +  omega x r,  omega = [0;0;1]
%
% eps and hence sigma are IDENTICAL to problemPatchTest: constant, and exactly
% representable. only the rotation differs, so any error here is purely the
% rotation coupling

    if nargin < 1 || isempty(lambda), lambda = 1; end
    if nargin < 2 || isempty(mu), mu = 1; end

    S0 = diag([2*mu + lambda, lambda, lambda]);

    problem.dim = 3;
    problem.name = sprintf('patch test + rigid rotation, lambda=%g mu=%g', lambda, mu);
    problem.mu     = @(x) mu * ones(1, size(x,2));
    problem.lambda = @(x) lambda * ones(1, size(x,2));
    problem.u      = @(x) [x(1,:) + 1 - x(2,:); x(1,:); zeros(1,size(x,2))];
    problem.sigma  = @(x) repmat(S0, 1, 1, size(x,2));
    problem.f      = @(x) zeros(3, size(x,2));

end
