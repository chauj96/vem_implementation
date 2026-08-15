function R_local = computeStressRecovery(cell_struct, face_struct, face_global_geom, B_local, P_local)
% construct the per-cell operator mapping global stress DOFs to the six
% projected polynomial stress coefficients
%
% P_E acts on local DOFs and x_local = T_E x_global, so
%
%   R_E = P_E T_E = P_E diag(s_E),   c_E = R_E sigma_h|_E

    faceLists = {cell_struct.faces};

    allFaces = cell2mat(cellfun(@(v) v(:), faceLists(:), 'UniformOutput', false));

    geomLists = cellfun(@(b) b.geom, B_local, 'UniformOutput', false);
    localGeom = [geomLists{:}];

    signs = computeFaceTransformations(allFaces, localGeom, face_global_geom(allFaces), face_struct);
    signs = signs(:);

    nCells = numel(cell_struct);
    R_local = cell(nCells,1);

    offset = 0;

    for e = 1:nCells

        m = 6 * numel(faceLists{e});

        R_local{e} = P_local{e} .* signs(offset + (1:m)).';

        offset = offset + m;

    end

end
