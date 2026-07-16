function [m20, m02, m11] = computeFaceMoments(X, xf, Qf)
% compute the second-order face moments

    nVertices = size(X, 1);

    % local coordinates of all face vertices
    Xshift = X' - xf;
    Xlocal = Qf' * Xshift;

    xtilde = Xlocal(1,:);
    ytilde = Xlocal(2,:);

    % planarity check
    normalCoordinates = Xlocal(3,:);

    faceScale = max(vecnorm(Xshift, 2, 1));

    if max(abs(normalCoordinates)) > 1e-10 * max(faceScale, eps)
        warning('Face vertices are not exactly coplanar.');
    end

    m20 = 0.0;
    m02 = 0.0;
    m11 = 0.0;

    for k = 1:nVertices

        kp1 = mod(k, nVertices) + 1;

        x1 = xtilde(k);
        y1 = ytilde(k);

        x2 = xtilde(kp1);
        y2 = ytilde(kp1);

        Jk = abs(x1*y2 - x2*y1);

        m20 = m20 + Jk * (x1^2 + x1*x2 + x2^2)/12.0;
        m02 = m02 + Jk * (y1^2 + y1*y2 + y2^2)/12.0;
        m11 = m11 + Jk * (x1*y1/12.0 + x1*y2/24.0 + x2*y1/24.0 + x2*y2/12.0);
    end
end