function signs = computeFaceTransformations(faces, localGeom, globalGeom, face_struct)
% construct the local-to-global stress DOF transformations with consistent face orientation
%
% T_f = diag(sign(n_L . n_G) * b),  b_j in {-1,+1} from  phi_j^L = b_j phi_j^G
%
%   faces      : P x 1 global face index of each element-face pair
%   localGeom  : 1 x P struct array of element-local face geometry
%   globalGeom : 1 x P struct array of global/reference face geometry
%   signs      : 6 x P, the diagonals of T_f
%
% batched over all element-face pairs at once

    % (component, quadrature point, face) values per vectorized chunk
    maxChunkValues = 2e6;

    faces = faces(:);
    P = numel(faces);

    signs = zeros(6, P);

    if P == 0
        return;
    end

    %% geometry as flat arrays

    t1L = reshape([localGeom.t1], 3, P);
    t2L = reshape([localGeom.t2], 3, P);
    nL = reshape([localGeom.n], 3, P);

    t1G = reshape([globalGeom.t1], 3, P);
    t2G = reshape([globalGeom.t2], 3, P);
    nG = reshape([globalGeom.n], 3, P);

    m20L = reshape([localGeom.m20], 1, P);
    m02L = reshape([localGeom.m02], 1, P);

    m20G = reshape([globalGeom.m20], 1, P);
    m02G = reshape([globalGeom.m02], 1, P);

    faceSub = face_struct(faces);

    Af = reshape([faceSub.area], 1, P);
    xf = reshape([faceSub.center], 3, P);

    quadPoints = {faceSub.quad_points};
    nQuad = cellfun(@(q) size(q,2), quadPoints);

    %% normal orientation

    normalDot = sum(nL .* nG, 1);

    bad = find(abs(abs(normalDot) - 1) >= 1e-10, 1);

    if ~isempty(bad)
        error('Local/global normals are not parallel on face %d.', faces(bad));
    end

    % local outward normal relative to global/reference normal
    normalSign = sign(normalDot);

    %% scaling of each basis function

    cL = [Af; Af; Af; ...
          Af .* sqrt(m20L + m02L); ...
          Af .* sqrt(m20L); ...
          Af .* sqrt(m02L)];

    cG = [Af; Af; Af; ...
          Af .* sqrt(m20G + m02G); ...
          Af .* sqrt(m20G); ...
          Af .* sqrt(m02G)];

    % max_q |phi_j^L -+ phi_j^G|
    errPlus = zeros(6, P);
    errMinus = zeros(6, P);

    %% phi_1..phi_3 constant over the face

    for j = 1:3

        switch j
            case 1, vL = t1L; vG = t1G;
            case 2, vL = t2L; vG = t2G;
            case 3, vL = nL;  vG = nG;
        end

        vL = vL ./ cL(j,:);
        vG = vG ./ cG(j,:);

        errPlus(j,:) = sqrt(sum((vL - vG).^2, 1));
        errMinus(j,:) = sqrt(sum((vL + vG).^2, 1));

    end

    %% phi_4..phi_6 linear over the face, grouped by quadrature point count

    for nq = unique(nQuad)

        group = find(nQuad == nq);

        % 3 nq values per face
        chunkSize = max(1, floor(maxChunkValues / (3*nq)));

        for c0 = 1:chunkSize:numel(group)

            sub = group(c0:min(c0+chunkSize-1, numel(group)));
            nSub = numel(sub);

            % x_q - x_f,  3 x nq x nSub
            qp = reshape([quadPoints{sub}], 3, nq, nSub);
            dx = qp - reshape(xf(:,sub), 3, 1, nSub);

            T1L = reshape(t1L(:,sub), 3, 1, nSub);
            T2L = reshape(t2L(:,sub), 3, 1, nSub);
            NL = reshape(nL(:,sub), 3, 1, nSub);

            T1G = reshape(t1G(:,sub), 3, 1, nSub);
            T2G = reshape(t2G(:,sub), 3, 1, nSub);
            NG = reshape(nG(:,sub), 3, 1, nSub);

            % xt = t1 . (x_q - x_f),  yt = t2 . (x_q - x_f)
            xtL = sum(T1L .* dx, 1);
            ytL = sum(T2L .* dx, 1);

            xtG = sum(T1G .* dx, 1);
            ytG = sum(T2G .* dx, 1);

            for j = 4:6

                sL = reshape(cL(j,sub), 1, 1, nSub);
                sG = reshape(cG(j,sub), 1, 1, nSub);

                switch j
                    case 4
                        PL = (T2L .* xtL - T1L .* ytL) ./ sL;
                        PG = (T2G .* xtG - T1G .* ytG) ./ sG;
                    case 5
                        PL = (NL .* xtL) ./ sL;
                        PG = (NG .* xtG) ./ sG;
                    case 6
                        PL = (NL .* ytL) ./ sL;
                        PG = (NG .* ytG) ./ sG;
                end

                % max over q of the euclidean norm
                errPlus(j,sub) = reshape(sqrt(max(sum((PL - PG).^2, 1), [], 2)), 1, nSub);
                errMinus(j,sub) = reshape(sqrt(max(sum((PL + PG).^2, 1), [], 2)), 1, nSub);

            end

        end

    end

    %% b_j

    scale = max(max(errPlus, errMinus), 1);
    tol = 1e-10 * scale;

    isPlus = errPlus < tol;                % phi_local =  phi_global
    isMinus = ~isPlus & errMinus < tol;    % phi_local = -phi_global

    matched = isPlus | isMinus;

    if ~all(matched, 'all')
        [j, p] = find(~matched, 1);
        error('Local/global basis cannot be matched by a sign change: face %d, basis %d.', faces(p), j);
    end

    basisSign = double(isPlus) - double(isMinus);

    signs = normalSign .* basisSign;

end
