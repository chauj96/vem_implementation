function problem = problemLinearStress(lambda, mu)
% u = [x^2/2; 0; 0]: linear stress, constant body force
%
% exercises the body-force path with a solution the discrete stress space can
% represent exactly on each face, isolating assembleBodyForce from the higher
% order behaviour of the 5.1 solution

    if nargin < 1 || isempty(lambda), lambda = 1; end
    if nargin < 2 || isempty(mu), mu = 1; end

    problem.dim = 3;
    problem.name = sprintf('linear stress (u = x^2/2), lambda=%g mu=%g', lambda, mu);
    problem.mu     = @(x) mu * ones(1, size(x,2));
    problem.lambda = @(x) lambda * ones(1, size(x,2));
    problem.u      = @(x) [0.5*x(1,:).^2; zeros(2, size(x,2))];

    % eps_xx = x, tr(eps) = x  ->  sigma = diag(2mu+lambda, lambda, lambda) * x
    problem.sigma  = @(x) reshape( ...
        [ (2*mu+lambda)*x(1,:); zeros(3,size(x,2)); lambda*x(1,:); zeros(3,size(x,2)); lambda*x(1,:) ], ...
        3, 3, []);

    % div(sigma) = [2mu+lambda; 0; 0]
    problem.f      = @(x) [(2*mu+lambda)*ones(1,size(x,2)); zeros(2,size(x,2))];

end
