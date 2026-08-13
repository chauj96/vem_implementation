function face_struct = assignFaceProperties(face_struct)
% find boundary faces and tagging boundary conditions

    tol = 1e-10;
    nFaces = numel(face_struct);
    
    for f = 1:nFaces
    
        face_struct(f).is_boundary = false;
        face_struct(f).is_dirichlet = false;
        face_struct(f).is_neumann = false;
    
        % exterior face
        if isscalar(face_struct(f).cells)
    
            face_struct(f).is_boundary = true;
    
            xf = face_struct(f).center(:);
    
            % Gamma_D: x = 0
            if abs(xf(1)) < tol
                face_struct(f).is_dirichlet = true;
    
            % Gamma_N: rest of boundary
            else
                face_struct(f).is_neumann = true;
            end
        end
    end

end