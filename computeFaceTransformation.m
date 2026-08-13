function Tf = computeFaceTransformation(f, localGeom, globalGeom, face_struct)
% construct the local-to-global stress DOF transformation with consistent face orientation

    xf = face_struct(f).center(:);
    Af = face_struct(f).area;
    qp = face_struct(f).quad_points;
    
    % local outward normal relative to global/reference normal
    normalSign = sign(dot(localGeom.n, globalGeom.n));
    
    assert(abs(abs(dot(localGeom.n, globalGeom.n)) - 1) < 1e-10, 'Local/global normals are not parallel on face %d.', f);
    
    Tf = zeros(6);
    
    for j = 1:6
    
        errPlus = 0;
        errMinus = 0;
    
        for q = 1:size(qp,2)
    
            xq = qp(:,q);
    
            % global/reference traction basis
            phiG = evaluateFaceBasis( ...
                j, xq, xf, globalGeom.Qf, ...
                globalGeom.t1, globalGeom.t2, globalGeom.n, ...
                Af, globalGeom.m20, globalGeom.m02, globalGeom.m11);
    
            % element-local traction basis
            phiL = evaluateFaceBasis( ...
                j, xq, xf, localGeom.Qf, ...
                localGeom.t1, localGeom.t2, localGeom.n, ...
                Af, localGeom.m20, localGeom.m02, localGeom.m11);
    
            errPlus = max(errPlus,  norm(phiL - phiG));
            errMinus = max(errMinus, norm(phiL + phiG));
        end
    
        scale = max([errPlus, errMinus, 1]);
        tol = 1e-10 * scale;
    
        % sign relating LOCAL BASIS to GLOBAL BASIS
        if errPlus < tol
            basisSign = 1;
    
        elseif errMinus < tol
            basisSign = -1;
    
        else
            error('Local/global basis cannot be matched by a sign change: face %d, basis %d.', f, j);
        end
    
        % stress DOF transformation also includes normal orientation
        Tf(j,j) = normalSign * basisSign;
    
    end

end