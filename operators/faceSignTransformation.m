function d = faceSignTransformation(f, localGeom, globalGeom)
% diagonal of the local -> global stress DOF transformation, in closed form
%
% Equivalent to diag(computeFaceTransformation(...)) but O(1): it derives the
% signs from the face frames instead of comparing the two bases at every
% quadrature point.
%
% The six face basis functions are
%
%   phi_1 = t1/Af                      phi_4 = (-yt t1 + xt t2)/(Af sqrt(m20+m02))
%   phi_2 = t2/Af                      phi_5 = (xt n)/(Af sqrt(m20))
%   phi_3 = n /Af                      phi_6 = (yt n)/(Af sqrt(m02))
%
% with xt = t1.(x - xf), yt = t2.(x - xf). If the local frame differs from the
% global one only by signs, t1_L = s1 t1_G, t2_L = s2 t2_G, n_L = s3 n_G, then
% xt_L = s1 xt_G and yt_L = s2 yt_G, so
%
%   basisSign = [s1, s2, s3, s1 s2, s1 s3, s2 s3]
%
% and m20, m02 are unchanged because they are even in xt, yt. The frames are
% verified to be sign-related; anything else is not a diagonal transformation
% and is rejected rather than silently mishandled.

    c = [localGeom.t1(:).' * globalGeom.t1(:), ...
         localGeom.t2(:).' * globalGeom.t2(:), ...
         localGeom.n(:).'  * globalGeom.n(:)];

    assert(all(abs(abs(c) - 1) < 1e-10), ...
        'Local/global face frames are not sign-related on face %d.', f);

    sg = sign(c);

    s1 = sg(1); s2 = sg(2); s3 = sg(3);

    % the in-plane moments must match, or the bases differ by more than a sign
    scale = max(globalGeom.m20 + globalGeom.m02, realmin);

    assert(abs(localGeom.m20 - globalGeom.m20) <= 1e-10 * scale && ...
           abs(localGeom.m02 - globalGeom.m02) <= 1e-10 * scale, ...
        'Local/global face moments differ on face %d.', f);

    basisSign = [s1, s2, s3, s1*s2, s1*s3, s2*s3];

    % normalSign in computeFaceTransformation is sign(n_L . n_G) = s3
    d = (s3 * basisSign).';

end
