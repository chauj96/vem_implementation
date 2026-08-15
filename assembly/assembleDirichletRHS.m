function rhs_D = assembleDirichletRHS(cell_struct, face_struct, face_global_geom)
% assemble Dirichlet boundary contribution to the stress RHS.
%
%   (rhs_D)_{f,j} = integral_f u_D * (tau_j n_out) dS
%
% stress DOFs are represented using the global/reference face convention.
% the orientation between the global face normal and the physical outward
% normal is corrected on each dirichlet boundary face.

    nFaces = numel(face_struct);
    nStressDofs = 6 * nFaces;
    
    rhs_D = zeros(nStressDofs,1);
       
    % for each boundary face, store its unique adjacent cell
    faceCell = zeros(nFaces,1);
    
    for e = 1:numel(cell_struct)
    
        elementFaces = cell_struct(e).faces;
    
        for lf = 1:numel(elementFaces)
    
            f = elementFaces(lf);
    
            % only boundary faces need an adjacent-cell lookup here
            if face_struct(f).is_dirichlet
                faceCell(f) = e;
            end
    
        end
    
    end
    
    % Dirichlet boundary contribution    
    dirichletFaces = find([face_struct.is_dirichlet]);
    
    for k = 1:numel(dirichletFaces)
    
        f = dirichletFaces(k);
    
        xf = face_struct(f).center(:);
        Af = face_struct(f).area;
    
        qp = face_struct(f).quad_points;
        w = face_struct(f).quad_weights;
    
        geom = face_global_geom(f);
    
        % physical outward orientation
        e = faceCell(f);
    
        assert(e ~= 0, 'Could not find adjacent cell for Dirichlet face %d.', f);
    
        xE = cell_struct(e).center(:);
    
        % direction from cell interior toward boundary face
        outwardDirection = xf - xE;
    
        % orientation of global/reference normal relative to the physical outward normal
        normalSign = sign(dot(geom.n, outwardDirection));
    
        assert(normalSign ~= 0, 'Could not determine outward orientation for face %d.', f);
 
        % Dirichlet integral: six RHS entries associated with this face
        faceRHS = zeros(6,1);
    
        for q = 1:numel(w)
    
            xq = qp(:,q);
    
            % prescribed displacement
            uD = [xq(1) + 1;
                  0;
                  0];
    
            % all six traction basis functions in the global/reference face convention
            PhiG = evaluateFaceBasis(xq, xf, geom.Qf, geom.t1, geom.t2, geom.n, Af, geom.m20, geom.m02);
    
            % convert global/reference orientation to physical outward orientation 
            PhiOut = normalSign * PhiG;
    
            % PhiOut' * uD gives all six integrands
            faceRHS = faceRHS + w(q) * (PhiOut.' * uD);
    
        end
        
        globalIdx = 6*(f-1) + (1:6);
    
        rhs_D(globalIdx) = faceRHS;
    
    end

end