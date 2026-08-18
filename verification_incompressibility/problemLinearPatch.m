function problem = problemLinearPatch(G, box, lambda, mu, name)
% linear-displacement patch test on an arbitrary mesh
%
%   u(x) = G (x - x0) / L
%
% with x0 the bounding-box centre and L its characteristic length, so the
% displacement is O(1) across the domain regardless of the mesh units.
%
% grad u = G/L is constant, hence
%
%   eps   = sym(G)/L                       constant
%   sigma = 2 mu eps + lambda tr(eps) I    constant
%   f     = div sigma = 0
%
% PASS CRITERION: sigma is constant and u is linear, so both lie in the discrete
% spaces exactly and every error norm is a consistency check, not an
% approximation error. the attainable floor is set by the linear solve, so the
% test passes at roughly one order of magnitude above the solver tolerance

    if nargin < 3 || isempty(lambda), lambda = 1; end
    if nargin < 4 || isempty(mu), mu = 1; end
    if nargin < 5 || isempty(name), name = 'linear patch'; end

    x0 = box.centre(:);
    L = box.length;

    Gs = G / L;
    eps0 = (Gs + Gs.') / 2;
    S0 = 2*mu*eps0 + lambda*trace(eps0)*eye(3);

    problem.dim = 3;
    problem.name = sprintf('%s, lambda=%g mu=%g', name, lambda, mu);
    problem.mu     = @(x) mu * ones(1, size(x,2));
    problem.lambda = @(x) lambda * ones(1, size(x,2));
    problem.u      = @(x) Gs * (x - x0);
    problem.sigma  = @(x) repmat(S0, 1, 1, size(x,2));
    problem.f      = @(x) zeros(3, size(x,2));

    problem.G = G;
    problem.sigmaConst = S0;

end
