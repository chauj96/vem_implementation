function problem = problemKeilegavlen51(dim, kappa, alpha)
% manufactured solution of Keilegavlen & Nordbotten (2017), section 5.1.1
%
%   dim   : 2 or 3. dim = 2 is the plane solution extended as u_z = 0 with no
%           z dependence, which is an exact 3D solution, run on a one-cell slab
%   kappa : heterogeneity contrast, chi = 1 in the upper corner (5.1.1)/(5.1.3)
%   alpha : lambda / mu
%
% with X = x-1/2 etc. and den = (1 - chi) + kappa*chi,
%
%   2D (5.1.2):  u = [ X^2 Y^2 ; -2/3 X Y^3 ; 0 ] / den
%   3D (5.1.4):  u = [ X^2Y^2Z^2 ; X^2Y^2Z^2 ; -2/3 (XY^2 + X^2Y) Z^3 ] / den
%
%   mu = den,  lambda = alpha*mu
%
% both solutions are divergence free, so tr(eps) = 0 and
%
%   sigma = 2 mu eps(u) = 2 eps(u_num)
%
% which is independent of the heterogeneity: sigma and the body force
% f = div(sigma) are smooth across the material interface.
%
% problem returns function handles, each taking a 3 x nq array of points:
%   u(x)      3 x nq displacement
%   sigma(x)  3 x 3 x nq stress
%   f(x)      3 x nq body force, f = div(sigma)
%   mu(x), lambda(x)   1 x nq

    if nargin < 2 || isempty(kappa), kappa = 1; end
    if nargin < 3 || isempty(alpha), alpha = 1; end

    problem.dim = dim;
    problem.kappa = kappa;
    problem.alpha = alpha;
    if kappa == 1
        % chi cancels: den == 1, so mu == 1 and lambda == alpha everywhere
        problem.name = sprintf('Keilegavlen 5.1 %dD, homogeneous, alpha=%g', dim, alpha);
    else
        problem.name = sprintf('Keilegavlen 5.1 %dD, kappa=%g, alpha=%g', dim, kappa, alpha);
    end

    if dim == 2
        problem.chi = @(x) double(min(x(1,:), x(2,:)) > 0.5);
    else
        problem.chi = @(x) double(min(min(x(1,:), x(2,:)), x(3,:)) > 0.5);
    end

    den = @(x) (1 - problem.chi(x)) + kappa * problem.chi(x);

    problem.mu = @(x) den(x);
    problem.lambda = @(x) alpha * den(x);

    if dim == 2
        problem.u     = @(x) u2d(x) ./ den(x);
        problem.sigma = @(x) sigma2d(x);
        problem.f     = @(x) f2d(x);
    else
        problem.u     = @(x) u3d(x) ./ den(x);
        problem.sigma = @(x) sigma3d(x);
        problem.f     = @(x) f3d(x);
    end

end

%% ------------------------------------------------------------------ 2D

function u = u2d(x)
    X = x(1,:) - 0.5;  Y = x(2,:) - 0.5;
    u = [ X.^2 .* Y.^2
         -(2/3) * X .* Y.^3
          zeros(size(X)) ];
end

function S = sigma2d(x)
% sigma = 2 eps(u_num);  eps_xx = 2XY^2, eps_yy = -2XY^2,
% eps_xy = X^2 Y - Y^3/3, all other components zero
    X = x(1,:) - 0.5;  Y = x(2,:) - 0.5;
    nq = numel(X);

    sxx =  4 * X .* Y.^2;
    syy = -4 * X .* Y.^2;
    sxy =  2 * (X.^2 .* Y - Y.^3 / 3);

    S = zeros(3,3,nq);
    S(1,1,:) = sxx;  S(2,2,:) = syy;
    S(1,2,:) = sxy;  S(2,1,:) = sxy;
end

function f = f2d(x)
% f = div(sigma)
    X = x(1,:) - 0.5;  Y = x(2,:) - 0.5;
    f = [ 2*(X.^2 + Y.^2)
         -4 * X .* Y
          zeros(size(X)) ];
end

%% ------------------------------------------------------------------ 3D

function u = u3d(x)
    X = x(1,:) - 0.5;  Y = x(2,:) - 0.5;  Z = x(3,:) - 0.5;
    u = [ X.^2 .* Y.^2 .* Z.^2
          X.^2 .* Y.^2 .* Z.^2
         -(2/3) * (X .* Y.^2 + X.^2 .* Y) .* Z.^3 ];
end

function S = sigma3d(x)
% sigma = 2 eps(u_num)
    X = x(1,:) - 0.5;  Y = x(2,:) - 0.5;  Z = x(3,:) - 0.5;
    nq = numel(X);

    exx =  2 * X .* Y.^2 .* Z.^2;
    eyy =  2 * X.^2 .* Y .* Z.^2;
    ezz = -2 * (X .* Y.^2 + X.^2 .* Y) .* Z.^2;

    exy = (X.^2 .* Y + X .* Y.^2) .* Z.^2;
    exz = X.^2 .* Y.^2 .* Z - (1/3) * (Y.^2 + 2*X.*Y) .* Z.^3;
    eyz = X.^2 .* Y.^2 .* Z - (1/3) * (2*X.*Y + X.^2) .* Z.^3;

    S = zeros(3,3,nq);
    S(1,1,:) = 2*exx;  S(2,2,:) = 2*eyy;  S(3,3,:) = 2*ezz;
    S(1,2,:) = 2*exy;  S(2,1,:) = 2*exy;
    S(1,3,:) = 2*exz;  S(3,1,:) = 2*exz;
    S(2,3,:) = 2*eyz;  S(3,2,:) = 2*eyz;
end

function f = f3d(x)
% f = div(sigma)
    X = x(1,:) - 0.5;  Y = x(2,:) - 0.5;  Z = x(3,:) - 0.5;

    common = 2 * (X.^2 .* Y.^2 + X.^2 .* Z.^2 + Y.^2 .* Z.^2);

    f = [ common
          common
         -4 * X .* Y .* (X + Y) .* Z - (4/3) * (X + Y) .* Z.^3 ];
end
