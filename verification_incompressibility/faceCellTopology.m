function topo = faceCellTopology(cell_struct, face_struct)
% face -> cell incidence graph, and the boundary it implies
%
% a facet with a single incident cell is external; one with two is internal.
% anything else means the mesh is not a valid watertight decomposition, so it
% is reported rather than silently accepted
%
%   topo.nCellsPerFace  nFaces x 1 incidence count
%   topo.isBoundary     nFaces x 1 logical, true where the count is 1
%   topo.faceCells      nFaces x 2, incident cell ids (0 where absent)
%   topo.orphan/.nonManifold   indices of faces with 0 or >2 incident cells

    nFaces = numel(face_struct);
    nCells = numel(cell_struct);

    % rebuild the incidence from the cell -> face lists rather than trusting
    % face_struct.cells, so the two representations are cross-checked
    faceCells = zeros(nFaces, 2);
    count = zeros(nFaces, 1);

    for e = 1:nCells
        for f = cell_struct(e).faces(:).'
            count(f) = count(f) + 1;
            if count(f) <= 2
                faceCells(f, count(f)) = e;
            end
        end
    end

    topo.nCellsPerFace = count;
    topo.isBoundary = (count == 1);
    topo.faceCells = faceCells;
    topo.orphan = find(count == 0);
    topo.nonManifold = find(count > 2);

    topo.nBoundary = sum(topo.isBoundary);
    topo.nInterior = sum(count == 2);

    % cross-check against what the mesh reader stored
    stored = cellfun(@numel, {face_struct.cells}).';
    topo.mismatch = find(stored ~= count);

end
