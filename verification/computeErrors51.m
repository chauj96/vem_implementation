function err = computeErrors51(cell_struct, face_struct, face_global_geom, sigma_h, u_h, problem)
% discrete L2 errors of Keilegavlen & Nordbotten (2017), (5.1.5) and (5.1.6)
%
%   e_D = [ sum_K m_K |u(x_K) - u_K|^2 ]^1/2 / [ sum_K m_K |u(x_K)|^2 ]^1/2
%   e_S = [ sum_s m_s |T(x_s) - T_s|^2 ]^1/2 / [ sum_s m_s |T(x_s)|^2 ]^1/2
%
% u_K is the discrete displacement at the cell barycenter: the rigid-motion
% rotation vanishes at x_E, so u_K is the translation part of the cell DOFs.
%
% the discrete stress unknown on a face is the traction, so the stress error
% is measured on tractions: with phi_1..phi_3 having face-mean t1, t2, n,
%
%   T_s = (1/|s|) integral_s sigma n dS = (c_1 t1 + c_2 t2 + c_3 n) / |s|
%
% compared against T(x_s) = sigma(x_s) n_s in the same reference orientation
%
% err.rotation is an extra diagnostic outside the paper's norms: the cell
% rotation DOF against the cell average of the exact rotation
%
%   omega_E = argmin_w  integral_E |u - alpha - w x (x - x_E)|^2
%
% which is exact for a rigid motion

    nCells = numel(cell_struct);
    nFaces = numel(face_struct);

    %% displacement, at cell barycenters
    xK = reshape([cell_struct.center], 3, nCells);
    mK = [cell_struct.volume].';

    uEx = problem.u(xK);                            % 3 x nCells
    uNum = reshape(u_h, 6, nCells);
    uNum = uNum(1:3, :);

    dU = sum((uEx - uNum).^2, 1).';
    nU = sum(uEx.^2, 1).';

    err.displacement = sqrt(mK.' * dU) / sqrt(mK.' * nU);

    %% traction, at face centers
    xS = reshape([face_struct.center], 3, nFaces);
    mS = [face_struct.area].';

    t1 = reshape([face_global_geom.t1], 3, nFaces);
    t2 = reshape([face_global_geom.t2], 3, nFaces);
    nn = reshape([face_global_geom.n],  3, nFaces);

    c = reshape(sigma_h, 6, nFaces);

    Tnum = (t1 .* c(1,:) + t2 .* c(2,:) + nn .* c(3,:)) ./ mS.';

    S = problem.sigma(xS);                          % 3 x 3 x nFaces
    Tex = squeeze(sum(S .* reshape(nn, 1, 3, nFaces), 2));

    dT = sum((Tex - Tnum).^2, 1).';
    nT = sum(Tex.^2, 1).';

    err.stress = sqrt(mS.' * dT) / sqrt(mS.' * nT);

    %% rotation, L2 projection of u onto RM(E)
    %
    % the discrete omega_E converges to the RM projection of u, not to the cell
    % average of curl(u)/2. those agree only when the cell inertia tensor is
    % isotropic (cubes); on tetrahedra or perturbed cells they differ at O(1),
    % so the projection is the correct reference. with R = [e1 x r, e2 x r, e3 x r]
    %
    %   (sum_q w_q R'R) omega = sum_q w_q R' u(x_q)
    omegaRef = zeros(3, nCells);

    for e = 1:nCells
        qp = cell_struct(e).quad_points;
        qw = cell_struct(e).quad_weights(:).';
        r = qp - cell_struct(e).center(:);
        uq = problem.u(qp);

        G = zeros(3,3);
        b = zeros(3,1);

        for q = 1:numel(qw)
            R = [      0,  r(3,q), -r(2,q);
                 -r(3,q),       0,  r(1,q);
                  r(2,q), -r(1,q),       0];
            G = G + qw(q) * (R.' * R);
            b = b + qw(q) * (R.' * uq(:,q));
        end

        omegaRef(:,e) = G \ b;
    end

    omegaNum = reshape(u_h, 6, nCells);
    omegaNum = omegaNum(4:6, :);

    dW = sum((omegaRef - omegaNum).^2, 1).';
    nW = sum(omegaRef.^2, 1).';

    refW = sqrt(mK.' * nW);

    % omegaRef comes from finite differences, so a curl-free field gives noise
    % rather than an exact zero. dividing by that noise would report a
    % meaningless blow-up, so fall back to the absolute error whenever the
    % reference rotation is negligible against the displacement scale
    uScale = sqrt(mK.' * nU);

    if refW > 1e-8 * max(uScale, realmin)
        err.rotation = sqrt(mK.' * dW) / refW;
        err.rotationIsRelative = true;
    else
        err.rotation = sqrt(mK.' * dW);
        err.rotationIsRelative = false;
    end

end

