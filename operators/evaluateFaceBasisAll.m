function Phi = evaluateFaceBasisAll(x, xf, Qf, t1, t2, n, Af, m20, m02)
% evaluate all six face traction basis functions at all quadrature points
%
%   x   : 3 x nq quadrature points
%   Phi : (3 nq) x 6, rows 3(q-1)+(1:3) hold Phi_q = evaluateFaceBasis(x(:,q), ...)

    nq = size(x, 2);

    t1 = t1(:);
    t2 = t2(:);
    n = n(:);

    % xt = t1 . (x - x_f),  yt = t2 . (x - x_f)
    local = Qf.' * (x - xf(:));

    xt = local(1,:);
    yt = local(2,:);

    % component rows of each quadrature point block
    i1 = 1:3:3*nq;
    i2 = 2:3:3*nq;
    i3 = 3:3:3*nq;

    Phi = zeros(3*nq, 6);

    % phi_1..phi_3 = t1/|f|, t2/|f|, n/|f|
    Phi(i1,1) = t1(1)/Af;  Phi(i2,1) = t1(2)/Af;  Phi(i3,1) = t1(3)/Af;
    Phi(i1,2) = t2(1)/Af;  Phi(i2,2) = t2(2)/Af;  Phi(i3,2) = t2(3)/Af;
    Phi(i1,3) = n(1)/Af;   Phi(i2,3) = n(2)/Af;   Phi(i3,3) = n(3)/Af;

    % phi_4 = (xt t2 - yt t1)/(|f| sqrt(m20+m02)),
    % phi_5 = xt n/(|f| sqrt(m20)),  phi_6 = yt n/(|f| sqrt(m02))
    c4 = (t2 .* xt - t1 .* yt) / (Af * sqrt(m20 + m02));
    c5 = (n .* xt) / (Af * sqrt(m20));
    c6 = (n .* yt) / (Af * sqrt(m02));

    Phi(:,4) = c4(:);
    Phi(:,5) = c5(:);
    Phi(:,6) = c6(:);

end
