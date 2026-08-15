function face_struct = computeFaceQuadrature(face_struct, V3)

  % precompute second-order triangle quadrature points and weights.

    alpha = 2/3;
    beta = 1/6;
    
    lambda = [ ...
        alpha beta  beta;
        beta  alpha beta;
        beta  beta  alpha];
    
    nFaces = numel(face_struct);

    for lf = 1:nFaces

        fC = face_struct(lf).center';  % barycenter of each face

        quad_points = [];
        quad_weights = [];
        
        face_verts = face_struct(lf).verts;

        for lv = 1:numel(face_verts)

            % choosing two consecutive vertices
            v1 = face_verts(lv);
            v2idx = mod(lv, numel(face_verts)) + 1;
            v2 = face_verts(v2idx);

            v1_coordinate = V3(v1, :);
            v2_coordinate = V3(v2, :);

            % define (decomposed) a single triangle
            triVerts = [ ...
                fC;
                v1_coordinate;
                v2_coordinate];

            % physical quadrature points
            Xq = lambda * triVerts;

            quad_points = [quad_points Xq'];

            % area of triangle
            edge1 = v1_coordinate - fC;
            edge2 = v2_coordinate - fC;

            areaTri = 0.5 * norm(cross(edge1, edge2));

            quad_weights = [quad_weights ...
                            areaTri / 3*ones(1,3)];
        end

        face_struct(lf).quad_points = quad_points;
        face_struct(lf).quad_weights = quad_weights;
        face_struct(lf).nQuad = size(quad_points, 2);
        face_struct(lf).nTri = face_struct(lf).nQuad / 3;
    end
end