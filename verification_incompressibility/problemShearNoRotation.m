function problem = problemShearNoRotation(mu, lambda)
% u = [x*y; x^2/2; 0]: off-diagonal stress but IDENTICALLY ZERO rotation
%
%   grad u is already symmetric, so curl(u) = 0
%   eps = [y, x, 0; x, 0, 0; 0, 0, 0],  tr(eps) = y
%   sigma = [ (2mu+lambda) y, 2 mu x, 0 ; 2 mu x, lambda y, 0 ; 0, 0, lambda y ]
%   f = div(sigma) = [0; 2 mu + lambda; 0]
%
% paired with problemShear (shear stress AND varying rotation) this separates
% the two candidate causes

    if nargin < 1 || isempty(mu), mu = 1; end
    if nargin < 2 || isempty(lambda), lambda = 1; end

    problem.dim = 3;
    problem.name = sprintf('shear stress, zero rotation, mu=%g lambda=%g', mu, lambda);
    problem.mu     = @(x) mu * ones(1, size(x,2));
    problem.lambda = @(x) lambda * ones(1, size(x,2));
    problem.u      = @(x) [x(1,:).*x(2,:); 0.5*x(1,:).^2; zeros(1,size(x,2))];

    problem.sigma  = @(x) reshape( [ ...
        (2*mu+lambda)*x(2,:);   2*mu*x(1,:);          zeros(1,size(x,2)); ...
        2*mu*x(1,:);            lambda*x(2,:);        zeros(1,size(x,2)); ...
        zeros(1,size(x,2));     zeros(1,size(x,2));   lambda*x(2,:) ], 3, 3, []);

    problem.f      = @(x) [zeros(1,size(x,2)); (2*mu+lambda)*ones(1,size(x,2)); zeros(1,size(x,2))];

end
