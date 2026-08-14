function Tf = computeFaceTransformation(f, localGeom, globalGeom, face_struct)
% construct the local-to-global stress DOF transformation with consistent face orientation

    xf = face_struct(f).center(:);
    Af = face_struct(f).area;
    qp = face_struct(f).quad_points;
    
    % normal orientation  
    normalDot = dot(localGeom.n, globalGeom.n);
    
    assert(abs(abs(normalDot) - 1) < 1e-10, 'Local/global normals are not parallel on face %d.', f);
    
    % local outward normal relative to global/reference normal
    normalSign = sign(normalDot);
    
    % compare local and global face bases
    errPlus = zeros(1,6);
    errMinus = zeros(1,6);
    
    for q = 1:size(qp,2)
    
        xq = qp(:,q);
    
        % global/reference face basis 
        PhiG = evaluateFaceBasis( ...
            xq, xf, ...
            globalGeom.Qf, ...
            globalGeom.t1, globalGeom.t2, globalGeom.n, ...
            Af, ...
            globalGeom.m20, globalGeom.m02);
    
        % element-local face basis 
        PhiL = evaluateFaceBasis( ...
            xq, xf, ...
            localGeom.Qf, ...
            localGeom.t1, localGeom.t2, localGeom.n, ...
            Af, ...
            localGeom.m20, localGeom.m02);

        % compare all six basis functions simultaneously
        plusError = sqrt(sum((PhiL - PhiG).^2, 1));
        minusError = sqrt(sum((PhiL + PhiG).^2, 1));
    
        errPlus = max(errPlus,  plusError);
        errMinus = max(errMinus, minusError);
    
    end
    
    % determine basis signs   
    scale = max([errPlus; errMinus; ones(1,6)], [], 1);
    tol = 1e-10 * scale;
    
    basisSign = zeros(1,6);
    
    for j = 1:6
    
        if errPlus(j) < tol(j)
    
            % phi_local = phi_global
            basisSign(j) = 1;
    
        elseif errMinus(j) < tol(j)
    
            % phi_local = -phi_global
            basisSign(j) = -1;
    
        else
    
            error('Local/global basis cannot be matched by a sign change: face %d, basis %d.', f, j);
    
        end
    
    end

    % local -> global stress DOF transformation
    Tf = diag(normalSign * basisSign);

end