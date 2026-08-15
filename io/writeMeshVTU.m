function writeMeshVTU(filename, V3, cell_struct, face_struct, cellData)
% Write a 3D polyhedral mesh to a VTU file for visualization in ParaView.
%
% cellData is an optional struct whose fields become named cell arrays, each
% nCells x nComp with nComp in {1, 3, 9} (scalar, vector, tensor).

    if nargin < 5 || isempty(cellData)
        cellData = struct();
    end

    nCells = numel(cell_struct);
    nPts = size(V3,1);

    VTK_POLYHEDRON = 42;

    connectivity = [];
    offsets = zeros(nCells,1);
    types = VTK_POLYHEDRON * ones(nCells,1,'uint8');

    % Polyhedron face connectivity
    faces = [];
    faceoffsets = zeros(nCells,1);

    off_conn = 0;
    off_face = 0;

    for c = 1:nCells

        fids = cell_struct(c).faces(:)';

        % Vertices used by this cell
        vids = [];

        for k = 1:numel(fids)
            vids = [vids, face_struct(fids(k)).verts(:)'];
        end

        vids = unique(vids, 'stable');   % MATLAB 1-based
        vids0 = vids - 1;                % VTK 0-based

        % Cell connectivity
        connectivity = [connectivity, vids0];

        off_conn = off_conn + numel(vids0);
        offsets(c) = off_conn;

        % VTK polyhedron face encoding:
        % [numFaces, nV(f1), v1, v2, ..., nV(f2), v1, v2, ..., ...]
        rec = numel(fids);

        for k = 1:numel(fids)

            v = face_struct(fids(k)).verts(:)' - 1;

            rec = [rec, numel(v), v];

        end

        faces = [faces, rec];

        off_face = off_face + numel(rec);
        faceoffsets(c) = off_face;

    end

    %% Write VTU file

    outputDir = fileparts(filename);

    if ~isempty(outputDir) && ~exist(outputDir, 'dir')
        mkdir(outputDir);
    end

    fid = fopen(filename, 'w');

    if fid == -1
        error('Could not open file: %s', filename);
    end

    fprintf(fid, '<?xml version="1.0"?>\n');
    fprintf(fid, '<VTKFile type="UnstructuredGrid" version="0.1" byte_order="LittleEndian">\n');

    fprintf(fid, '<UnstructuredGrid>\n');
    fprintf(fid, '<Piece NumberOfPoints="%d" NumberOfCells="%d">\n', nPts, nCells);

    %% Cell data

    writeCellData(fid, cellData, nCells);

    %% Points

    fprintf(fid, '<Points>\n');

    fprintf(fid, '<DataArray type="Float64" NumberOfComponents="3" format="ascii">\n');

    fprintf(fid, '%.15g %.15g %.15g\n', V3');

    fprintf(fid, '</DataArray>\n');
    fprintf(fid, '</Points>\n');

    %% Cells

    fprintf(fid, '<Cells>\n');

    % Cell-to-vertex connectivity
    fprintf(fid, '<DataArray type="Int32" Name="connectivity" format="ascii">\n');

    fprintf(fid, '%d ', connectivity);
    fprintf(fid, '\n</DataArray>\n');

    % Connectivity offsets
    fprintf(fid, '<DataArray type="Int32" Name="offsets" format="ascii">\n');

    fprintf(fid, '%d ', offsets);
    fprintf(fid, '\n</DataArray>\n');

    % Cell types
    fprintf(fid, '<DataArray type="UInt8" Name="types" format="ascii">\n');

    fprintf(fid, '%d ', types);
    fprintf(fid, '\n</DataArray>\n');

    % Polyhedron faces
    fprintf(fid, '<DataArray type="Int32" Name="faces" format="ascii">\n');

    fprintf(fid, '%d ', faces);
    fprintf(fid, '\n</DataArray>\n');

    % Polyhedron face offsets
    fprintf(fid, '<DataArray type="Int32" Name="faceoffsets" format="ascii">\n');

    fprintf(fid, '%d ', faceoffsets);
    fprintf(fid, '\n</DataArray>\n');

    fprintf(fid, '</Cells>\n');

    %% Finish file

    fprintf(fid, '</Piece>\n');
    fprintf(fid, '</UnstructuredGrid>\n');
    fprintf(fid, '</VTKFile>\n');

    fclose(fid);

    names = fieldnames(cellData);

    if isempty(names)
        fprintf('\n[# Mesh written to %s]\n', filename);
    else
        fprintf('\n[# Mesh written to %s, cell fields: %s]\n', filename, strjoin(names.', ', '));
    end

end

%% helper
function writeCellData(fid, cellData, nCells)
% write every field of cellData as a named VTK cell array

    names = fieldnames(cellData);

    if isempty(names)
        return;
    end

    % the first field of each kind becomes the active one in ParaView
    active = struct('Scalars', '', 'Vectors', '', 'Tensors', '');

    nComponents = zeros(numel(names),1);

    for k = 1:numel(names)

        A = cellData.(names{k});

        if size(A,1) ~= nCells
            error('Cell field "%s" has %d rows, expected %d.', names{k}, size(A,1), nCells);
        end

        nComponents(k) = size(A,2);

        switch nComponents(k)
            case 1, kind = 'Scalars';
            case 3, kind = 'Vectors';
            case 9, kind = 'Tensors';
            otherwise
                error('Cell field "%s" has %d components, expected 1, 3 or 9.', names{k}, nComponents(k));
        end

        if isempty(active.(kind))
            active.(kind) = names{k};
        end

    end

    attributes = '';

    for kind = {'Scalars','Vectors','Tensors'}

        if ~isempty(active.(kind{1}))
            attributes = [attributes, sprintf(' %s="%s"', kind{1}, active.(kind{1}))]; %#ok<AGROW>
        end

    end

    fprintf(fid, '<CellData%s>\n', attributes);

    for k = 1:numel(names)

        A = cellData.(names{k});

        fprintf(fid, '<DataArray type="Float64" Name="%s" NumberOfComponents="%d" format="ascii">\n', ...
            names{k}, nComponents(k));

        % one cell per line
        format = [repmat('%.15g ', 1, nComponents(k)), '\n'];

        fprintf(fid, format, A.');

        fprintf(fid, '</DataArray>\n');

    end

    fprintf(fid, '</CellData>\n');

end