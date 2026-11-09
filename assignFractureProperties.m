function [cell_struct, face_struct, fracture] = assignFractureProperties(cell_struct, face_struct, fractureInfo)
% Defines the fracture geometry, identifies its faces and adjacent cells, and establishes a consistent orientation for the interface formulation.

    nFaces = numel(face_struct);
    nCells = numel(cell_struct);

    %% Initialize fracture geometry

    fracture.faces = fractureInfo.faces(:)';
    fracture.nFaces = numel(fracture.faces);

    % Normal points from minus to plus
    fracture.normal = fractureInfo.normal(:);
    fracture.tangent = fractureInfo.tangent(:);

    fracture.minus_cells = fractureInfo.minusCells(:);
    fracture.plus_cells = fractureInfo.plusCells(:);

    fracture.psi = fractureInfo.psi;
    fracture.b = fractureInfo.b;
    fracture.area = fractureInfo.area;

    %% Initialize face properties

    for f = 1:nFaces
        face_struct(f).is_fracture = false;
        face_struct(f).fracture_id = 0;
    end

    %% Validate and assign fracture faces

    assert(fracture.nFaces > 0, 'No fracture faces identified.');

    assert(numel(unique(fracture.faces)) == fracture.nFaces, 'Duplicate fracture face IDs.');

    assert(all(fracture.faces >= 1 & fracture.faces <= nFaces), 'Invalid fracture face IDs.');

    assert(numel(fracture.minus_cells) == fracture.nFaces && numel(fracture.plus_cells) == fracture.nFaces, 'Fracture cell count mismatch.');

    assert(abs(norm(fracture.normal) - 1) < 1e-10, 'Fracture normal must be a unit vector.');

    for k = 1:fracture.nFaces

        f = fracture.faces(k);

        cm = fracture.minus_cells(k);
        cp = fracture.plus_cells(k);

        assert(cm >= 1 && cm <= nCells && cp >= 1 && cp <= nCells && cm ~= cp, 'Invalid fracture cells at face %d.', f);

        neighbors = face_struct(f).cells(:);

        % Fracture must be an internal face
        assert(numel(neighbors) == 2, 'Fracture face %d is not internal.', f);

        % Check adjacent cells
        assert(isequal(sort(neighbors), sort([cm; cp])), 'Fracture neighbor mismatch at face %d.', f);

        % Check minus/plus orientation
        xf = face_struct(f).center(:);
        xm = cell_struct(cm).center(:);
        xp = cell_struct(cp).center(:);

        assert(dot(xm - xf, fracture.normal) < 0 && dot(xp - xf, fracture.normal) > 0, 'Incorrect fracture orientation at face %d.', f);

        face_struct(f).is_fracture = true;
        face_struct(f).fracture_id = k;

    end

    %% Summary

    fprintf('\n[# FRACTURE GEOMETRY]\n');
    fprintf('  Inclination: %.2f degrees\n', fracture.psi);
    fprintf('  Fracture length: %.4f\n', 2*fracture.b);
    fprintf('  Fracture faces: %d\n', fracture.nFaces);
    fprintf('  Fracture area: %.6f\n', fracture.area);

end
