function [G, fractureInfo] = unitCubeFracture(Nx, Ny, Nz, xFrac, yRange, zRange, meshType)
% generate a conforming unit cube mesh with a planar fracture.
%
% fracture:
%   x = xFrac
%   y in yRange
%   z in zRange

    Lx = 1;
    Ly = 1;
    Lz = 1;

    tol = 1e-10;

    % generate bulk mesh
    G = cartGrid([Nx Ny Nz], [Lx Ly Lz], 'meshType', meshType);
    G = computeGeometry(G);

    % initialize fracture face tags
    G.faces.is_fracture = false(G.faces.num, 1);

    fractureFaces = [];
    minusCells = [];
    plusCells = [];

    % global fracture orientation: minus -> plus
    normal = [1; 0; 0];

    for f = 1:G.faces.num

        neighbors = G.faces.neighbors(f, :);

        % fracture must be an internal face
        if any(neighbors == 0)
            continue;
        end

        xf = G.faces.centroids(f, :)';

        % check fracture plane
        if abs(xf(1) - xFrac) > tol
            continue;
        end

        % check finite fracture region
        if xf(2) < yRange(1)-tol || xf(2) > yRange(2)+tol
            continue;
        end

        if xf(3) < zRange(1)-tol || xf(3) > zRange(2)+tol
            continue;
        end

        % check normal direction
        nf = G.faces.normals(f, :)';
        nf = nf / norm(nf);

        if abs(abs(dot(nf, normal)) - 1) > tol
            continue;
        end

        % check fracture face vertices
        nodeIds = G.faces.nodes(G.faces.nodePos(f):G.faces.nodePos(f+1)-1);
        X = G.nodes.coords(nodeIds, :);

        if any(abs(X(:,1) - xFrac) > tol)
            error('Fracture face is not coplanar.');
        end

        if any(X(:,2) < yRange(1)-tol | X(:,2) > yRange(2)+tol)
            error('Fracture boundary is not aligned with the mesh.');
        end

        if any(X(:,3) < zRange(1)-tol | X(:,3) > zRange(2)+tol)
            error('Fracture boundary is not aligned with the mesh.');
        end

        % identify minus and plus cells
        c1 = neighbors(1);
        c2 = neighbors(2);

        x1 = G.cells.centroids(c1, :)';
        x2 = G.cells.centroids(c2, :)';

        d1 = dot(x1 - xf, normal);
        d2 = dot(x2 - xf, normal);

        if d1 < 0 && d2 > 0

            cm = c1;
            cp = c2;

        elseif d2 < 0 && d1 > 0

            cm = c2;
            cp = c1;

        else

            error('Incorrect fracture orientation at face %d.', f);

        end

        % store fracture face information
        G.faces.is_fracture(f) = true;

        fractureFaces(end+1, 1) = f;
        minusCells(end+1, 1) = cm;
        plusCells(end+1, 1) = cp;

    end

    nF = numel(fractureFaces);

    if nF == 0
        error('No fracture faces identified.');
    end

    % store fracture information
    fractureInfo = struct();

    fractureInfo.faces = fractureFaces;
    fractureInfo.nFaces = nF;

    fractureInfo.minusCells = minusCells;
    fractureInfo.plusCells = plusCells;

    fractureInfo.cells = [minusCells, plusCells];

    fractureInfo.normal = normal;

    fractureInfo.center = [xFrac; mean(yRange); mean(zRange)];

    fractureInfo.xFrac = xFrac;
    fractureInfo.yRange = yRange;
    fractureInfo.zRange = zRange;

    fractureInfo.area = sum(G.faces.areas(fractureFaces));

    % check total fracture area
    expectedArea = diff(yRange) * diff(zRange);

    if abs(fractureInfo.area - expectedArea) > tol * max(expectedArea, 1)
        error('Fracture area mismatch. Check mesh alignment.');
    end

    fprintf('\n[# FRACTURE GEOMETRY]\n');
    fprintf('  Fracture plane: x = %.3f\n', xFrac);
    fprintf('  Fracture faces: %d\n', nF);
    fprintf('  Fracture area: %.6f\n', fractureInfo.area);

end