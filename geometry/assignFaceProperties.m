function face_struct = assignFaceProperties(face_struct, Lx, Ly, Lz)
% Identify exterior boundary faces and classify them by their geometric location (xmin, xmax, ymin, ymax, zmin, zmax).
% Initialize boundary condition flags without imposing any BCs. Internal faces, including fracture faces, are excluded.

    tol = 1e-7 * max([Lx, Ly, Lz, 1]);

    nFaces = numel(face_struct);

    for f = 1:nFaces

        face_struct(f).is_boundary = false;
        face_struct(f).is_displacement_bc = false;
        face_struct(f).is_traction_bc = false;

        face_struct(f).boundary_side = '';

        if numel(face_struct(f).cells) ~= 1
            continue;
        end

        face_struct(f).is_boundary = true;

        x = face_struct(f).center(:);

        if abs(x(1) + Lx/2) < tol
            side = 'xmin';
        elseif abs(x(1) - Lx/2) < tol
            side = 'xmax';
        elseif abs(x(2) + Ly/2) < tol
            side = 'ymin';
        elseif abs(x(2) - Ly/2) < tol
            side = 'ymax';
        elseif abs(x(3)) < tol
            side = 'zmin';
        elseif abs(x(3) - Lz) < tol
            side = 'zmax';
        else
            error('Unrecognized exterior face %d.', f);
        end

        face_struct(f).boundary_side = side;
    end

    fprintf('\n[# BOUNDARY FACE PROPERTIES]\n');
    fprintf('  Exterior boundary faces: %d\n', sum([face_struct.is_boundary]));

end

