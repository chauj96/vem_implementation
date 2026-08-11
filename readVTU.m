function [cell_struct, face_struct, V3, cells3D] = readVTU(filename)
% read VTU mesh files and return geometry information

    doc = xmlread(filename);
    
    % vertices
    pts = doc.getElementsByTagName('Points');
    da = pts.item(0).getElementsByTagName('DataArray').item(0);
    
    txt = char(da.getFirstChild.getData);
    V3 = sscanf(txt,'%f');
    V3 = reshape(V3,3,[])';
    
    % necessary arrays
    cellCenter = readDataArray(doc,'cellCenter');
    cellVolume = readDataArray(doc,'cellVolume');
    
    faceCenter = readDataArray(doc,'faceCenter');
    faceNormal = readDataArray(doc,'faceNormal');
    faceArea = readDataArray(doc,'faceArea');
    
    faceCells = readDataArray(doc,'faceCells','int64');
    
    cellFacesFlat = readDataArray(doc,'cellFaces_flat','int64');
    cellFaceOffsets = readDataArray(doc,'cellFaceOffsets','int64');
    
    faceVertsFlat = readDataArray(doc,'faceVerts_flat','int64');
    faceVertOffsets = readDataArray(doc,'faceVertOffsets','int64');
    
    nFaces = size(faceCenter,1);
    nCells = size(cellCenter,1);
    
    face_struct = struct([]);
    cell_struct = struct([]);
    cells3D = cell(nCells,1);
    
    %% face_struct  
    for f = 1:nFaces
    
        if f == 1
            ids = 1:faceVertOffsets(1);
        else
            ids = faceVertOffsets(f-1)+1:faceVertOffsets(f);
        end
    
        verts = faceVertsFlat(ids);
    
        cells = faceCells(f,:);
        cells = cells(cells>=0) + 1; % indexing!
    
        face_struct(f).cells = cells(:)';
        face_struct(f).verts = verts(:)+1;
        face_struct(f).center = faceCenter(f,:)';
        face_struct(f).normal = faceNormal(f,:)';
        face_struct(f).area = faceArea(f);
    
    end
    
    %% cell_struct  
    for c = 1:nCells
    
        if c == 1
            ids = 1:cellFaceOffsets(1);
        else
            ids = cellFaceOffsets(c-1)+1:cellFaceOffsets(c);
        end
    
        faces = cellFacesFlat(ids) + 1;
        xc = cellCenter(c,:)';
        signs = zeros(numel(faces),1);
        normals = zeros(numel(faces),3);
    
        verts_local = [];
    
        for k = 1:numel(faces)
    
            f = faces(k);
    
            xf = face_struct(f).center;
            nf = face_struct(f).normal;
    
            if dot(nf, xf - xc) < 0
                signs(k) = -1;
                normals(k,:) = (-nf)';
            else
                signs(k) = 1;
                normals(k,:) = nf';
            end
    
            verts_local = [verts_local; face_struct(f).verts];
    
        end
    
        cell_struct(c).faces = faces(:)';
        cell_struct(c).center = xc;
        cell_struct(c).volume = cellVolume(c);
    
        cell_struct(c).faces_orientation = signs';
        cell_struct(c).face_normals = normals;
    
        cells3D{c} = unique(verts_local,'stable')';
    
    end
    
    fprintf('\n');
    fprintf('3D mesh model info:\n');
    fprintf(' %d vertices\n',size(V3,1));
    fprintf(' %d cells\n',nCells);
    fprintf(' %d faces\n',nFaces);
    
    checkVolumeConsistency(cell_struct,face_struct,nCells);

end

%% helper 1
function checkVolumeConsistency(cell_struct, face_struct, nCells)
% verify cell geometry via the MFD consistency condition.
%
%   Computes the tensor T = sum_f [ (s * n_f * A_f) * (x_f - x_c)^T ].
%   For a valid polyhedron, T must equal V_c * I (where I is 3x3 identity).

    inexact = 0;
    max_rel_err = 0;
    fail_list = [];

    fprintf('\n=== MRST Grid Convert:: MFD Volume & Tensor Consistency Check ===\n');

    for c = 1:nCells
        cell_faces = cell_struct(c).faces;
        face_signs = cell_struct(c).faces_orientation;
        Vc = cell_struct(c).volume;
        xc = cell_struct(c).center(:);

        T = zeros(3,3);

        for k = 1:numel(cell_faces)
            f = cell_faces(k);
            s  = face_signs(k); % +1 or -1 outward sign
            xf = face_struct(f).center(:);
            nf = face_struct(f).normal(:); % unit normal
            af = face_struct(f).area;

            n_out = s * nf; % Outward area-weighted normal (3x1)
            T = T + af * n_out * (xf - xc)';
        end

        % Trace/div-theorem volume
        V_div = trace(T) / 3;

        % Check MFD tensor equation
        rel_err = norm(T - Vc * eye(3), 'fro') / max(abs(Vc), eps);
        rel_vol_err = abs(V_div - Vc) / max(abs(Vc), eps);
        max_rel_err = max(max_rel_err, rel_err);

        if rel_err > 1e-10
            inexact = inexact + 1;
            fail_list(end+1,:) = [c, rel_err, Vc, V_div, rel_vol_err]; 
        end
    end

    fprintf('  Cells checked           : %d\n', nCells);
    fprintf('  Max relative tensor err : %.3e\n', max_rel_err);

    if inexact == 0
        fprintf('  All cells PASS (tol 1e-10).\n');
    else
        fprintf('  FAILED cells (rel_err > 1e-10): %d / %d\n', inexact, nCells);

        vol_err = fail_list(:,5);
        fprintf('\n  Volume discrepancy (|V_div - V_stored|/V_stored) in failed cells:\n');
        fprintf('    min/mean/max: %.3e / %.3e / %.3e\n', ...
                min(vol_err), mean(vol_err), max(vol_err));

        [~, ord] = sort(fail_list(:,2), 'descend');
        fprintf('\n  Worst 5 cells (ci, rel_T_err, V_stored, V_div):\n');
        for k = 1:min(5, size(fail_list,1))
            r = fail_list(ord(k),:);
            fprintf('    cell %4d  rel_err=%.3e  V_stored=%.4e  V_div=%.4e\n', ...
                    r(1), r(2), r(3), r(4));
        end
    end
    fprintf('-----------------------------------------------------\n\n');
end

%% helper 2
function data = readDataArray(doc,name,dtype)
% extracts the specified XML data array, converts its contents to either
% double or int64 and reshapes multicomponent arrays according to the
% NumberOfComponents attribute.

    if nargin < 3
        dtype = 'double';
    end

    arrays = doc.getElementsByTagName('DataArray');

    for k = 0:arrays.getLength-1

        node = arrays.item(k);

        if node.hasAttribute('Name')

            if strcmp(char(node.getAttribute('Name')),name)

                txt = char(node.getFirstChild.getData);

                nComp = 1;

                if node.hasAttribute('NumberOfComponents')
                    nComp = str2double(char(node.getAttribute('NumberOfComponents')));
                end

                switch lower(dtype)

                    case 'double'
                        data = sscanf(txt,'%f');

                    case 'int64'
                        data = int64(sscanf(txt,'%ld'));

                    otherwise
                        error('Unsupported data type.');

                end

                if nComp > 1
                    data = reshape(data,nComp,[])';
                end

                return

            end
        end
    end

    error('Could not find DataArray "%s".',name)

end