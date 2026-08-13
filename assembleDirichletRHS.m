function rhs_D = assembleDirichletRHS(face_struct, face_global_geom)
% assemble Dirichlet boundary contribution to the stress RHS:
%
%   (rhs_D)_{f,j} = integral_f u_D · phi^G_{f,j} dS
%
% where phi^G uses the global/reference face convention.

    nFaces = numel(face_struct);
    nStressDofs = 6 * nFaces;
    
    rhs_D = zeros(nStressDofs, 1);
    
    for f = 1:nFaces
    
        % only Dirichlet boundary faces
        if ~face_struct(f).is_dirichlet
            continue;
        end
    
        xf = face_struct(f).center(:);
        Af = face_struct(f).area;
    
        qp = face_struct(f).quad_points;
        w = face_struct(f).quad_weights;
    
        geom = face_global_geom(f);
    
        % six traction basis functions on face f
        for j = 1:6
    
            val = 0.0;
    
            for q = 1:numel(w)
    
                xq = qp(:,q);
    
                % prescribed displacement (given)
                uD = [xq(1) + 1;
                      0;
                      0];
    
                % global/reference traction basis
                phiG = evaluateFaceBasis(j, xq, xf, geom.Qf, geom.t1, geom.t2, geom.n, Af, geom.m20, geom.m02, geom.m11);
                val = val + w(q) * dot(uD, phiG);
            end
    
            % global stress DOF corresponding to face f, basis j
            globalIdx = 6*(f-1) + j;
    
            rhs_D(globalIdx) = val;
        end
    end

end