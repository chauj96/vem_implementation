function writeMeshVTU(filename, V3, cell_struct, face_struct)
% Write a 3D polyhedral mesh to a VTU file for visualization in ParaView.

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

    fprintf('\n[# Mesh written to %s]\n', filename);

end