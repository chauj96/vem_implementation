function [cell_struct, face_struct, V3, cells3D] = simplexMeshBox(n, L, perturb, seed)
% tetrahedral mesh of the box [0,L(1)] x [0,L(2)] x [0,L(3)]
%
% each hexahedron is split into six tetrahedra sharing the main diagonal
% (Kuhn decomposition), matching the paper's "simplex grids created by
% dividing squares". triangular faces are planar under any perturbation.
%
% arguments as in cartesianMeshBox

    if nargin < 2 || isempty(L), L = [1 1 1]; end
    if nargin < 3 || isempty(perturb), perturb = 0; end
    if nargin < 4 || isempty(seed), seed = 0; end

    n = n(:).'; if isscalar(n), n = [n n n]; end
    nx = n(1); ny = n(2); nz = n(3);

    [X, Y, Z] = ndgrid(linspace(0,L(1),nx+1), linspace(0,L(2),ny+1), linspace(0,L(3),nz+1));
    V3 = [X(:), Y(:), Z(:)];
    V3 = perturbVertices(V3, n, L, perturb, seed);

    vid = @(i,j,k) i + (nx+1)*(j-1) + (nx+1)*(ny+1)*(k-1);

    % Kuhn decomposition of the unit cube, corners indexed (di,dj,dk)
    corner = @(i,j,k,d) vid(i+d(1), j+d(2), k+d(3));
    paths = [0 0 0; 1 0 0; 1 1 0; 1 1 1
             0 0 0; 1 1 0; 0 1 0; 1 1 1
             0 0 0; 0 1 0; 0 1 1; 1 1 1
             0 0 0; 0 1 1; 0 0 1; 1 1 1
             0 0 0; 0 0 1; 1 0 1; 1 1 1
             0 0 0; 1 0 1; 1 0 0; 1 1 1];

    nTet = 6*nx*ny*nz;
    tets = zeros(nTet,4);
    t = 0;
    for k = 1:nz, for j = 1:ny, for i = 1:nx
        for s = 1:6
            t = t + 1;
            d = paths(4*(s-1)+(1:4), :);
            tets(t,:) = [corner(i,j,k,d(1,:)), corner(i,j,k,d(2,:)), ...
                         corner(i,j,k,d(3,:)), corner(i,j,k,d(4,:))];
        end
    end, end, end

    % triangular faces of every tet, identified by their sorted vertex triple
    localTri = [1 2 3; 1 2 4; 1 3 4; 2 3 4];
    triAll = zeros(4*nTet, 3);
    for s = 1:4
        triAll((s-1)*nTet + (1:nTet), :) = tets(:, localTri(s,:));
    end

    [~, ia, ic] = unique(sort(triAll, 2), 'rows', 'stable');

    nFaces = numel(ia);
    faceVerts = num2cell(triAll(ia,:), 2);
    faceCells = zeros(nFaces,2);

    tetOfEntry = repmat((1:nTet).', 4, 1);
    slot = ones(nFaces,1);
    for e = 1:numel(ic)
        f = ic(e);
        if slot(f) <= 2
            faceCells(f, slot(f)) = tetOfEntry(e);
            slot(f) = slot(f) + 1;
        end
    end

    cellFaces = cell(nTet,1);
    for f = 1:nFaces
        for c = faceCells(f, faceCells(f,:) > 0)
            cellFaces{c}(end+1) = f; %#ok<AGROW>
        end
    end

    [cell_struct, face_struct, cells3D] = buildMeshStructs(V3, faceVerts, faceCells, cellFaces);

end
