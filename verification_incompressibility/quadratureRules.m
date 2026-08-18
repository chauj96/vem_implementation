function [bary, w] = quadratureRules(shape, degree)
% barycentric quadrature rules of arbitrary degree on the reference triangle
% or tetrahedron, built from Gauss-Legendre points via the Duffy transform
%
%   shape  : 'tri' or 'tet'
%   degree : polynomial degree to integrate exactly
%
%   bary : nq x (3 or 4) barycentric coordinates
%   w    : nq x 1 weights summing to 1, so physical weights are w * measure
%
% the Duffy map collapses the square/cube onto the simplex; its Jacobian is
% polynomial, so a Gauss rule of n = ceil((degree+3)/2) points per direction
% integrates the transformed degree-`degree` integrand exactly

    n = max(1, ceil((degree + 3) / 2));
    [x, gw] = gaussLegendre01(n);

    switch lower(shape)

        case 'tri'
            [U, V] = ndgrid(x, x);
            [WU, WV] = ndgrid(gw, gw);
            U = U(:); V = V(:); W = WU(:) .* WV(:);

            a = U;
            b = V .* (1 - U);
            jac = (1 - U);

            w = W .* jac;
            w = w / sum(w);
            bary = [1 - a - b, a, b];

        case 'tet'
            [U, V, T] = ndgrid(x, x, x);
            [WU, WV, WT] = ndgrid(gw, gw, gw);
            U = U(:); V = V(:); T = T(:);
            W = WU(:) .* WV(:) .* WT(:);

            a = U;
            b = V .* (1 - U);
            c = T .* (1 - U) .* (1 - V);
            jac = (1 - U).^2 .* (1 - V);

            w = W .* jac;
            w = w / sum(w);
            bary = [1 - a - b - c, a, b, c];

        otherwise
            error('Unknown shape "%s".', shape);
    end

end

function [x, w] = gaussLegendre01(n)
% n-point Gauss-Legendre nodes and weights on [0,1] (Golub-Welsch)

    if n == 1
        x = 0.5; w = 1;
        return;
    end

    k = 1:n-1;
    beta = k ./ sqrt(4*k.^2 - 1);
    J = diag(beta, 1) + diag(beta, -1);

    [V, D] = eig(J);
    [xr, idx] = sort(diag(D));
    wr = 2 * (V(1,idx).^2).';

    x = (xr + 1) / 2;      % map [-1,1] -> [0,1]
    w = wr / 2;

end
