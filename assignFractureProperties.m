function [cell_struct, face_struct, fracture] = assignFractureProperties(cell_struct, face_struct, xFrac, yRange, zRange)
% identify fracture faces and adjacent bulk cells.
%
% geometry-only function, independent of elasticity/flow.
%
% fracture:
%   x = xFrac
%   y in yRange
%   z in zRange
%
% example:
%   assignFractureProperties(..., 0.5, [0.25 0.75], [0.25 0.75])

    if nargin < 4 || isempty(yRange)
        yRange = [0.25, 0.75];
    end

    if nargin < 5 || isempty(zRange)
        zRange = [0.25, 0.75];
    end

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

        xf = face_struct(f).center(:);

        % check fracture plane
        if abs(xf(1) - xFrac) > tol
            continue;
        end

        % check finite fracture region
        if xf(2) < yRange(1)-tol || xf(2) > yRange(2)+tol || xf(3) < zRange(1)-tol || xf(3) > zRange(2)+tol
            continue;
        end

        % check normal direction
        nf = face_struct(f).normal(:);
        nf = nf / norm(nf);

        if abs(abs(nf(1)) - 1) > tol
            continue;
        end

        face_struct(f).is_fracture = true;
        fracture.faces(end+1) = f;
    end

    fracture.nFaces = numel(fracture.faces);

    % global orientation: minus -> plus
    fracture.normal = [1; 0; 0];

    fracture.minus_cells = zeros(fracture.nFaces, 1);
    fracture.plus_cells = zeros(fracture.nFaces, 1);

    for k = 1:fracture.nFaces

        f = fracture.faces(k);
        neighbors = face_struct(f).cells;

        c1 = neighbors(1);
        c2 = neighbors(2);

        if cell_struct(c1).center(1) < cell_struct(c2).center(1)

            fracture.minus_cells(k) = c1;
            fracture.plus_cells(k) = c2;
        else
            fracture.minus_cells(k) = c2;
            fracture.plus_cells(k) = c1;
        end

        face_struct(f).fracture_id = k;
    end

    fracture.xFrac = xFrac;
    fracture.yRange = yRange;
    fracture.zRange = zRange;

    fprintf('\n[# FRACTURE GEOMETRY]\n');
    fprintf('  Fracture plane: x = %.3f\n', xFrac);
    fprintf('  Fracture faces: %d\n', fracture.nFaces);
end
