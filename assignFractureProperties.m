function [cell_struct, face_struct, fracture] = assignFractureProperties(cell_struct, face_struct, xFrac)

    tol = 1e-10;
    nFaces = numel(face_struct);

    fracture.faces = [];

    for f = 1:nFaces

        face_struct(f).is_fracture = false;
        face_struct(f).fracture_id = 0;

        % fracture must be an internal face
        if numel(face_struct(f).cells) ~= 2
            continue;
        end

        % check whether face lies on x = xFrac
        if abs(face_struct(f).center(1) - xFrac) > tol
            continue;
        end

        % require normal parallel to x-axis
        nf = face_struct(f).normal(:);

        if abs(abs(nf(1)) - 1) > tol
            continue;
        end

        face_struct(f).is_fracture = true;

        fracture.faces(end+1) = f;
    end

    fracture.nFaces = numel(fracture.faces);

    % consistent orientation: normal points from minus to plus
    fracture.normal = [1; 0; 0];

    fracture.minus_cells = zeros(fracture.nFaces, 1);
    fracture.plus_cells = zeros(fracture.nFaces, 1);

    for k = 1:fracture.nFaces

        f = fracture.faces(k);
        neighbors = face_struct(f).cells;

        c1 = neighbors(1);
        c2 = neighbors(2);

        x1 = cell_struct(c1).center(1);
        x2 = cell_struct(c2).center(1);

        if x1 < x2
            fracture.minus_cells(k) = c1;
            fracture.plus_cells(k) = c2;
        else
            fracture.minus_cells(k) = c2;
            fracture.plus_cells(k) = c1;
        end

        face_struct(f).fracture_id = k;
    end

    fprintf('\n[# FRACTURE GEOMETRY]\n');
    fprintf('  Fracture plane: x = %.3f\n', xFrac);
    fprintf('  Fracture faces: %d\n', fracture.nFaces);
end