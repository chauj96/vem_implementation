function [cell_struct, face_struct, cells3D] = buildMeshStructs(V3, faceVerts, faceCells, cellFaces)
% build cell_struct / face_struct / cells3D from raw connectivity
%
%   V3        : nVerts x 3 vertex coordinates
%   faceVerts : nFaces x 1 cell, ordered vertex loop of each face
%   faceCells : nFaces x 2, adjacent cell ids (0 where there is no neighbour)
%   cellFaces : nCells x 1 cell, face ids of each cell
%
% all geometry is computed from V3, so the same connectivity can be reused
% after the vertices are perturbed. faces are fan-triangulated from their
% vertex mean, cells are split into tetrahedra on those triangles

    nFaces = numel(faceVerts);
    nCells = numel(cellFaces);

    face_struct = struct('cells', cell(nFaces,1), 'verts', [], ...
                         'center', [], 'normal', [], 'area', []);

    faceCentroid = zeros(nFaces,3);
    faceNormal = zeros(nFaces,3);
    faceArea = zeros(nFaces,1);

    for f = 1:nFaces

        p = V3(faceVerts{f}, :);
        m = size(p,1);

        c0 = mean(p, 1);

        % fan triangles (c0, p_k, p_{k+1}) : area vectors and centroids
        kp1 = [2:m, 1];
        a = p - c0;
        b = p(kp1,:) - c0;

        cr = 0.5 * cross(a, b, 2);            % triangle area vectors
        triArea = sqrt(sum(cr.^2, 2));
        triCent = (c0 + p + p(kp1,:)) / 3;

        areaVec = sum(cr, 1);
        A = norm(areaVec);

        faceArea(f) = A;
        faceNormal(f,:) = areaVec / A;
        faceCentroid(f,:) = sum(triArea .* triCent, 1) / sum(triArea);

        face_struct(f).verts = faceVerts{f}(:);
        face_struct(f).center = faceCentroid(f,:).';
        face_struct(f).normal = faceNormal(f,:).';
        face_struct(f).area = A;

        nb = faceCells(f, faceCells(f,:) > 0);
        face_struct(f).cells = nb(:).';

    end

    cell_struct = struct('center', cell(nCells,1), 'volume', [], 'faces', [], ...
                         'faces_orientation', [], 'face_normals', []);
    cells3D = cell(nCells,1);

    for c = 1:nCells

        fl = cellFaces{c}(:).';
        nf = numel(fl);

        % interior reference point
        x0 = mean(faceCentroid(fl,:), 1);

        % outward orientation of each face relative to this cell
        s = sign(sum(faceNormal(fl,:) .* (faceCentroid(fl,:) - x0), 2));
        s(s == 0) = 1;

        vol = 0;
        momentSum = zeros(1,3);
        vids = [];

        for k = 1:nf

            f = fl(k);
            p = V3(faceVerts{f}, :);
            m = size(p,1);
            c0 = mean(p, 1);
            kp1 = [2:m, 1];

            % tetrahedra (x0, c0, p_k, p_{k+1}), signed then oriented outward
            A = c0 - x0;
            B = p - x0;
            C = p(kp1,:) - x0;

            vt = s(k) * sum(cross(repmat(A, m, 1), B, 2) .* C, 2) / 6;

            ct = (x0 + c0 + p + p(kp1,:)) / 4;

            vol = vol + sum(vt);
            momentSum = momentSum + sum(vt .* ct, 1);

            vids = [vids; faceVerts{f}(:)]; %#ok<AGROW>

        end

        cell_struct(c).volume = vol;
        cell_struct(c).center = (momentSum / vol).';
        cell_struct(c).faces = fl;
        cell_struct(c).faces_orientation = s(:).';
        cell_struct(c).face_normals = faceNormal(fl,:) .* s;

        cells3D{c} = unique(vids, 'stable').';

    end

end
