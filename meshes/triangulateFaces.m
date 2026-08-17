function [faceVerts, faceCells, cellFaces] = triangulateFaces(faceVerts, faceCells, cellFaces)
% split every face with more than three vertices into a triangle fan
%
% a perturbed hexahedral mesh has non-planar quadrilateral faces, which the
% face frame and face moments do not model. splitting each face into triangles
% makes every facet planar by construction; cells become general polyhedra,
% which the VEM handles directly (a perturbed hex becomes a 12-face polyhedron)

    nFacesOld = numel(faceVerts);

    newVerts = {};
    newCells = [];
    childOf = cell(nFacesOld,1);

    for f = 1:nFacesOld

        v = faceVerts{f}(:).';
        m = numel(v);

        if m <= 3
            newVerts{end+1} = v;               %#ok<AGROW>
            newCells(end+1,:) = faceCells(f,:); %#ok<AGROW>
            childOf{f} = numel(newVerts);
            continue;
        end

        % fan from the first vertex
        idx = zeros(1, m-2);
        for k = 2:m-1
            newVerts{end+1} = [v(1), v(k), v(k+1)];  %#ok<AGROW>
            newCells(end+1,:) = faceCells(f,:);      %#ok<AGROW>
            idx(k-1) = numel(newVerts);
        end
        childOf{f} = idx;

    end

    faceVerts = newVerts(:);
    faceCells = newCells;

    for c = 1:numel(cellFaces)
        cellFaces{c} = [childOf{cellFaces{c}}];
    end

end
