function face_global_geom = computeGlobalFaceGeometry(face_struct, V3)
% construct global/reference geometry for each face

    nFaces = numel(face_struct);
    face_global_geom = struct([]);
    
    for f = 1:nFaces
    
        xf = face_struct(f).center(:);
    
        vertexIds = face_struct(f).verts(:);
        X = V3(vertexIds, :);
    
        % global/reference face normal
        n = face_struct(f).normal(:);
        n = n / norm(n);
    
        % global/reference face frame
        [t1, t2, Qf] = constructFaceFrame(X, n);
    
        % global/reference face moments
        [m20, m02, m11] = computeFaceMoments(X, xf, Qf);
    
        % save the global quantities
        face_global_geom(f).t1 = t1;
        face_global_geom(f).t2 = t2;
        face_global_geom(f).n = n;
        face_global_geom(f).Qf = Qf;
        face_global_geom(f).m20 = m20;
        face_global_geom(f).m02 = m02;
        face_global_geom(f).m11 = m11;
    
    end

end