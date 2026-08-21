function G = projectionBasisBatch(x, xE)
% evaluateProjectionBasis at every quadrature point at once
%
% x is 3-by-nq. G is (3 nq)-by-6, rows grouped by point: rows 3(q-1)+(1:3) hold
% the 3-by-6 block of point q. Batching removes one function call and one small
% matrix build per quadrature point.

    r = x - xE;

    xr = r(1,:).';
    yr = r(2,:).';
    zr = r(3,:).';

    nq = size(r,2);
    z = zeros(nq,1);

    s = 1/sqrt(2);

    G = zeros(3*nq, 6);

    % shear ordering matches projToMatrices: 4 -> yz, 5 -> xz, 6 -> xy
    G(1:3:end,:) = [xr,  z,   z,   z,     s*zr, s*yr];
    G(2:3:end,:) = [z,   yr,  z,   s*zr,  z,    s*xr];
    G(3:3:end,:) = [z,   z,   zr,  s*yr,  s*xr, z   ];

end
