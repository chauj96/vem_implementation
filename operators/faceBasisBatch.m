function Phi = faceBasisBatch(x, xf, Qf, t1, t2, n, Af, m20, m02)
% evaluateFaceBasis at every quadrature point at once
%
% x is 3-by-nq. Phi is (3 nq)-by-6, rows grouped by point: rows 3(q-1)+(1:3)
% hold the 3-by-6 block of point q.

    local = Qf.' * (x - xf);

    xt = local(1,:);
    yt = local(2,:);

    nq = size(x,2);

    t1 = t1(:); t2 = t2(:); n = n(:);

    Phi = [ ...
        reshape(repmat(t1, 1, nq), [], 1) / Af, ...
        reshape(repmat(t2, 1, nq), [], 1) / Af, ...
        reshape(repmat(n,  1, nq), [], 1) / Af, ...
        reshape(t2*xt - t1*yt, [], 1) / (Af*sqrt(m20+m02)), ...
        reshape(n*xt, [], 1) / (Af*sqrt(m20)), ...
        reshape(n*yt, [], 1) / (Af*sqrt(m02))];

end
