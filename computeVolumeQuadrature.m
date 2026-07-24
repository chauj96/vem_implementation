function cell_struct = computeVolumeQuadrature(cell_struct, face_struct, V3)

    % precompute second-order tetrahedral quadrature points and weights for
    % every polyhedral element.
    
    alpha = 0.5854101966249685;
    beta  = 0.1381966011250153;
    
    lambda = [ ...
        alpha beta  beta  beta;
        beta  alpha beta  beta;
        beta  beta  alpha beta;
        beta  beta  beta  alpha];
    
    nCells = numel(cell_struct);
    
    for c = 1:nCells
    
        xE = cell_struct(c).center(:);
    
        quad_points  = [];
        quad_weights = [];
    
        cell_faces = cell_struct(c).faces;
    
        for lf = 1:numel(cell_faces)
    
            f = cell_faces(lf);
    
            xf = face_struct(f).center(:);
    
            verts = face_struct(f).verts(:);
            nVerts = numel(verts);
    
            % triangulate face
            for k = 1:nVerts
    
                kp1 = mod(k, nVerts) + 1;
                vi   = V3(verts(k), :)';
                vip1 = V3(verts(kp1), :)';
    
                % tetrahedron vertices
                tetVerts = [ ...
                    xE';
                    xf';
                    vi';
                    vip1'];
    
                % physical quadrature points
                Xq = lambda * tetVerts;
    
                quad_points = [quad_points Xq'];
    
                % tetrahedron volume
                J = [xf - xE, vi - xE, vip1 - xE];
    
                volumeTet = abs(det(J)) / 6;
    
                quad_weights = [quad_weights ...
                                volumeTet / 4*ones(1,4)];
            end
    
        end
    
        cell_struct(c).quad_points = quad_points;
        cell_struct(c).quad_weights = quad_weights;
        cell_struct(c).nQuad = size(quad_points, 2);
        cell_struct(c).nTet = cell_struct(c).nQuad / 4;
    
    end

end