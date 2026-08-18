function problem = problemShear(mu, lambda)
% u = [0; x^2/2; 0]: same polynomial degree as problemLinearStress, but with a
% NON-ZERO rotation
%
%   eps_xy = x/2, tr(eps) = 0  ->  sigma_xy = sigma_yx = mu x
%   curl(u)/2 = [0; 0; x/2]                  (problemLinearStress has curl u = 0)
%   f = div(sigma) = [0; mu; 0]
%
% controlled comparison: if the stress converges here, the rotation coupling is
% sound and the 5.1 failure is about polynomial degree; if it fails, the
% rotation coupling is the culprit

    if nargin < 1 || isempty(mu), mu = 1; end
    if nargin < 2 || isempty(lambda), lambda = 1; end

    problem.dim = 3;
    problem.name = sprintf('shear (u = [0; x^2/2; 0]), mu=%g lambda=%g', mu, lambda);
    problem.mu     = @(x) mu * ones(1, size(x,2));
    problem.lambda = @(x) lambda * ones(1, size(x,2));
    problem.u      = @(x) [zeros(1,size(x,2)); 0.5*x(1,:).^2; zeros(1,size(x,2))];

    problem.sigma  = @(x) reshape( ...
        [ zeros(1,size(x,2)); mu*x(1,:); zeros(1,size(x,2)); ...
          mu*x(1,:); zeros(5,size(x,2)) ], 3, 3, []);

    problem.f      = @(x) [zeros(1,size(x,2)); mu*ones(1,size(x,2)); zeros(1,size(x,2))];

end
