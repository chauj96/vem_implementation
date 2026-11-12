function fracture = buildFractureSurfaceMesh(fracture, face_struct, V3)
% build a 2D fracture mesh from existing 3D fracture faces.
%
% each fracture face is a 2D cell, each unique fracture edge supports one tangential flux DOF.
%
% Assumptions:
%   - conforming planar fracture mesh
%   - face vertices are ordered along the polygon boundary
%   - fracture faces are convex (Cartesian benchmark)
%
% edge orientation:
%   global edge co-normal is defined using its first adjacent
%   fracture face. local signs convert this orientation into
%   outward co-normal orientation for each fracture face.

    nF = fracture.nFaces;

    if nF == 0
        error('No fracture faces found.');
    end

    %% 1. Collect fracture edges

    allEdges = zeros(0, 2);
    faceEdgeRows = cell(nF, 1);

    fracture.face_centers = zeros(nF, 3);
    fracture.face_areas = zeros(nF, 1);

    for k = 1:nF

        f = fracture.faces(k);
        verts = face_struct(f).verts(:);
        nv = numel(verts);

        % Consecutive vertices form polygon edges
        pairs = [verts, verts([2:nv, 1])];

        % Unique edge identification ignores vertex ordering
        pairs = sort(pairs, 2);

        firstRow = size(allEdges, 1) + 1;
        faceEdgeRows{k} = firstRow:(firstRow + nv - 1);

        allEdges = [allEdges; pairs];

        fracture.face_centers(k,:) = face_struct(f).center(:)';

        fracture.face_areas(k) = face_struct(f).area;
    end

    %% 2. Unique edges and connectivity

    [fracture.edges, ~, edgeMap] = unique(allEdges, 'rows', 'stable');

    nE = size(fracture.edges, 1);
    fracture.nEdges = nE;

    fracture.face_edges = cell(nF, 1);
    fracture.edge_faces = cell(nE, 1);

    for k = 1:nF

        eids = edgeMap(faceEdgeRows{k});
        fracture.face_edges{k} = eids(:)';

        for j = 1:numel(eids)
            e = eids(j);
            fracture.edge_faces{e}(end+1) = k;
        end
    end

    %% 3. Edge geometry

    fracture.edge_centers = zeros(nE, 3);
    fracture.edge_lengths = zeros(nE, 1);

    for e = 1:nE

        v1 = V3(fracture.edges(e,1), :);
        v2 = V3(fracture.edges(e,2), :);

        fracture.edge_centers(e,:) = (v1 + v2)/2;
        fracture.edge_lengths(e) = norm(v2 - v1);

        if fracture.edge_lengths(e) <= 0
            error('Degenerate fracture edge %d.', e);
        end
    end

    %% 4. Boundary and internal edges

    nAdjacent = cellfun(@numel, fracture.edge_faces);

    if any(nAdjacent > 2)
        error('Non-manifold fracture edge detected.');
    end

    fracture.boundary_edges = find(nAdjacent == 1);
    fracture.internal_edges = find(nAdjacent == 2);

    % Fracture perimeter boundary conditions
    fracture.is_boundary_edge = false(nE,1);
    fracture.is_boundary_edge(fracture.boundary_edges) = true;
    
    fracture.is_flux_bc_edge = false(nE,1);
    fracture.is_flux_bc_edge(fracture.boundary_edges) = true;
    
    fracture.is_pressure_bc_edge = false(nE,1);

    %% 5. Global co-normal and local orientation signs

    % global fracture plane normal
    nGamma = fracture.normal(:);
    nGamma = nGamma / norm(nGamma);

    fracture.edge_normals = zeros(nE, 3);
    fracture.face_edge_signs = cell(nF, 1);

    for k = 1:nF
        fracture.face_edge_signs{k} = zeros(1, numel(fracture.face_edges{k}));
    end

    for e = 1:nE

        % tangent vector along the edge
        v1 = V3(fracture.edges(e,1), :)';
        v2 = V3(fracture.edges(e,2), :)';

        t = (v2 - v1) / fracture.edge_lengths(e);

        % candidate co-normal in fracture plane
        nu = cross(nGamma, t);
        nu = nu / norm(nu);

        % fix global orientation using first adjacent face
        k0 = fracture.edge_faces{e}(1);

        xf = fracture.face_centers(k0,:)';
        xe = fracture.edge_centers(e,:)';

        d = xe - xf;

        %pProject onto fracture plane
        d = d - dot(d, nGamma)*nGamma;

        if dot(nu, d) < 0
            nu = -nu;
        end

        fracture.edge_normals(e,:) = nu';

        % determine local outward orientation for each face
        adjacent = fracture.edge_faces{e};

        for a = 1:numel(adjacent)

            k = adjacent(a);
            xf = fracture.face_centers(k,:)';

            d = xe - xf;
            d = d - dot(d, nGamma)*nGamma;

            s = sign(dot(nu, d));

            if s == 0
                error('Ambiguous edge orientation.');
            end

            loc = find(fracture.face_edges{k} == e);

            fracture.face_edge_signs{k}(loc) = s;
        end
    end

    %% 6. Check orientation consistency

    for e = fracture.internal_edges(:)'

        adjacent = fracture.edge_faces{e};

        k1 = adjacent(1);
        k2 = adjacent(2);

        j1 = find(fracture.face_edges{k1} == e);
        j2 = find(fracture.face_edges{k2} == e);

        s1 = fracture.face_edge_signs{k1}(j1);
        s2 = fracture.face_edge_signs{k2}(j2);

        if s1 ~= -s2
            error('Inconsistent orientation at edge %d.', e);
        end
    end

    fprintf('\n[# FRACTURE SURFACE MESH]\n');
    fprintf('  Fracture faces: %d\n', nF);
    fprintf('  Unique edges:   %d\n', nE);
    fprintf('  Internal edges: %d\n', numel(fracture.internal_edges));
    fprintf('  Boundary edges: %d\n', numel(fracture.boundary_edges));
end
