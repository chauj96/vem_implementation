function [G, fractureInfo] = inclinedFracture(Lx, Ly, Lz, h, Nz, psi, b)
% generate a 2D fracture-conforming PEBI mesh and extrude to 3D.
%
% domain:
%   [-Lx/2,Lx/2] x [-Ly/2,Ly/2] x [0,Lz]
%
% fracture:
%   centered at (0,0)
%   length = 2*b
%   angle = psi degrees
%
% outputs:
%   G                  - MRST 3D grid
%   fractureInfo.faces - 3D fracture face IDs
%   fractureInfo.cells - Neighboring cell pairs
%   fractureInfo.plusCells / minusCells - Oriented neighbors

    %% 1. Input checks and fracture geometry

    assert(Lx > 0 && Ly > 0 && Lz > 0 && h > 0 && b > 0);
    assert(Nz >= 1 && Nz == round(Nz));

    theta = deg2rad(psi);

    tangent = [cos(theta); sin(theta); 0];
    normal  = [-sin(theta); cos(theta); 0];

    p1 = [-b*cos(theta), -b*sin(theta)];
    p2 = [ b*cos(theta),  b*sin(theta)];

    assert(all(all(abs([p1; p2]) < [Lx/2, Ly/2])), 'Fracture endpoints must lie inside the domain.');

    %% 2. Generate 2D PEBI mesh

    fractureLine = {[p1; p2] + [Lx/2, Ly/2]};

    G2D = pebiGrid2D(h, [Lx, Ly], 'faceConstraints', fractureLine, 'interpolateFC', true);

    % center the domain at the origin
    G2D.nodes.coords(:,1) = G2D.nodes.coords(:,1) - Lx/2;
    G2D.nodes.coords(:,2) = G2D.nodes.coords(:,2) - Ly/2;

    G2D = computeGeometry(G2D);

    assert(G2D.cells.num > 0, '2D grid is empty.');

    %% 3. Identify fracture edges in 2D using MRST tags

    assert(isfield(G2D.faces, 'tag'), 'MRST did not generate face tags.');

    faces2D = find(G2D.faces.tag ~= 0);

    assert(~isempty(faces2D), 'No tagged fracture edges in 2D mesh.');

    neighbors2D = G2D.faces.neighbors(faces2D,:);

    assert(all(neighbors2D(:) > 0), 'Some tagged fracture edges are boundary faces.');

    fractureLength2D = sum(G2D.faces.areas(faces2D));

    lengthTol = 1e-4 * max(2*b,1);

    assert(abs(fractureLength2D - 2*b) < lengthTol, 'Tagged fracture length does not match 2*b.');

    %% 4. Extrude 2D mesh into 3D

    G = makeLayeredGrid(G2D, Nz);

    assert(G.cells.num == G2D.cells.num * Nz, '3D extrusion failed.');

    % Scale z coordinates to [0,Lz]
    z = G.nodes.coords(:,3);

    zmin = min(z);
    zmax = max(z);

    assert(zmax > zmin, 'Extrusion has zero thickness.');

    G.nodes.coords(:,3) = Lz * (z-zmin)/(zmax-zmin);

    G = computeGeometry(G);

    %% 5. Match 2D fracture edges to 3D faces
    %
    % a 3D fracture face must:
    %  - be an internal face
    %  - have vertices lying on a tagged 2D fracture edge
    %  - span a nonzero distance in z

    tol = 1e-7 * max([Lx,Ly,Lz,1]);

    fractureFaces = [];

    for f = 1:G.faces.num

        % only internal faces
        if any(G.faces.neighbors(f,:) == 0)
            continue;
        end

        nodeRange = G.faces.nodePos(f):G.faces.nodePos(f+1)-1;
        nodeIds = G.faces.nodes(nodeRange);

        xyz = G.nodes.coords(nodeIds,:);

        % skip horizontal faces
        if max(xyz(:,3))-min(xyz(:,3)) < tol
            continue;
        end

        xy = xyz(:,1:2);

        % check whether all XY vertices lie on one tagged 2D edge
        isFractureFace = false;

        for k = 1:numel(faces2D)

            f2 = faces2D(k);

            edgeRange = G2D.faces.nodePos(f2): G2D.faces.nodePos(f2+1)-1;

            edgeNodes = G2D.faces.nodes(edgeRange);
            edgeXY = G2D.nodes.coords(edgeNodes,:);

            % each 3D face vertex must match a 2D edge endpoint
            D = pdist2(xy, edgeXY);

            if all(min(D,[],2) < tol)
                isFractureFace = true;
                break;
            end
        end

        if isFractureFace
            fractureFaces(end+1,1) = f;
        end
    end

    expectedFaces = numel(faces2D)*Nz;

    assert(numel(fractureFaces) == expectedFaces, 'Expected %d fracture faces, found %d.', expectedFaces, numel(fractureFaces));

    %% 6. Identify neighboring cells and orient plus/minus

    fractureCells = G.faces.neighbors(fractureFaces,:);

    plusCells = zeros(numel(fractureFaces),1);
    minusCells = zeros(numel(fractureFaces),1);

    for k = 1:numel(fractureFaces)

        c1 = fractureCells(k,1);
        c2 = fractureCells(k,2);

        x1 = G.cells.centroids(c1,:);
        x2 = G.cells.centroids(c2,:);

        % Determine which cell is on positive normal side
        d1 = dot(x1, normal');
        d2 = dot(x2, normal');

        if d1 > d2
            plusCells(k) = c1;
            minusCells(k) = c2;
        else
            plusCells(k) = c2;
            minusCells(k) = c1;
        end
    end

    %% 7. Store fracture information

    fractureInfo.center = [0;0;0];
    fractureInfo.tangent = tangent;
    fractureInfo.normal = normal;

    fractureInfo.psi = psi;
    fractureInfo.b = b;
    fractureInfo.endpoints = [p1;p2];

    fractureInfo.faces2D = faces2D;
    fractureInfo.faces = fractureFaces;
    fractureInfo.nFaces = numel(fractureFaces);

    fractureInfo.cells = fractureCells;
    fractureInfo.plusCells = plusCells;
    fractureInfo.minusCells = minusCells;

    fractureInfo.area = sum(G.faces.areas(fractureFaces));

    expectedArea = 2*b*Lz;

    assert(abs(fractureInfo.area-expectedArea) < 1e-4*max(expectedArea,1), 'Fracture area does not match 2*b*Lz.');

    %% 

    fprintf('\n[# INCLINED FRACTURE MESH]\n');
    fprintf('  Inclination: %.2f degrees\n', psi);
    fprintf('  Fracture length: %.4f\n', 2*b);
    fprintf('  2D cells: %d\n', G2D.cells.num);
    fprintf('  3D cells: %d\n', G.cells.num);
    fprintf('  3D faces: %d\n', G.faces.num);
    fprintf('  Fracture faces (2D): %d\n', numel(faces2D));
    fprintf('  Fracture faces (3D): %d\n', fractureInfo.nFaces);
    fprintf('  Fracture area: %.6f\n', fractureInfo.area);
    fprintf('  Z thickness: %.4f\n', Lz);

end





