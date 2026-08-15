function [K_global, B_global, K_cons_global, K_stab_global] = assembleGlobalMatrices(cell_struct, face_struct, face_global_geom, B_local, K_cons, K_stab)
% assemble global mixed-VEM matrices using a common global face convention
%
%   K = sum_E  T_E' K_E T_E,     B = sum_E  B_E T_E,     K = K^cons + K^stab
%
% T_E = diag(s_E), s_E in {-1,+1}^{6 n_F(E)}, hence
%
%   (T_E' K_E T_E)_ij = s_i s_j (K_E)_ij,     (B_E T_E)_ij = (B_E)_ij s_j

    nCells = numel(cell_struct);
    nFaces = numel(face_struct);

    nStressDofs = 6 * nFaces;
    nDispDofs = 6 * nCells;

    % entries per vectorized chunk
    maxChunkEntries = 4e6;

    %% element sizes and local -> global stress DOF map

    faceLists = {cell_struct.faces};

    nElementFaces = cellfun(@numel, faceLists);
    nElementFaces = nElementFaces(:);

    % m_E = 6 n_F(E)
    blockSize = 6 * nElementFaces;

    nLocalDofs = sum(blockSize);
    nKEntries = sum(blockSize.^2);
    nBEntries = 6 * nLocalDofs;

    % offset of cell E in the concatenated local numbering
    cellOffset = [0; cumsum(blockSize)];

    % (E,lf) -> f
    allFaces = cell2mat(cellfun(@(v) v(:), faceLists(:), 'UniformOutput', false));

    B_matrix = cellfun(@(b) b.matrix, B_local, 'UniformOutput', false);

    % local DOF 6(lf-1)+i  ->  global DOF 6(f-1)+i
    stressDofMap = reshape(6*(allFaces(:).' - 1) + (1:6).', [], 1);

    %% s_E, the diagonal of T_E

    % ordered as allFaces
    geomLists = cellfun(@(b) b.geom, B_local, 'UniformOutput', false);
    localGeom = [geomLists{:}];

    assert(isequal([localGeom.face].', allFaces), ...
        'Element face lists and local face geometry are out of sync.');

    signs = computeFaceTransformations(allFaces, localGeom, face_global_geom(allFaces), face_struct);

    assert(all(abs(signs) == 1, 'all'), 'Invalid face transformation.');

    signs = signs(:);

    %% triplet construction, grouped by n_F(E)

    splitParts = nargout > 2;

    % K^cons and K^stab share (IK,JK)
    IK = zeros(nKEntries,1);
    JK = zeros(nKEntries,1);
    VKc = zeros(nKEntries,1);

    if splitParts
        VKs = zeros(nKEntries,1);
    end

    IB = zeros(nBEntries,1);
    JB = zeros(nBEntries,1);
    VB = zeros(nBEntries,1);

    kPtr = 0;
    bPtr = 0;

    for nf = unique(nElementFaces).'

        group = find(nElementFaces == nf);

        m = 6 * nf;

        % m^2 entries per cell
        chunkSize = max(1, floor(maxChunkEntries / m^2));

        for c0 = 1:chunkSize:numel(group)

            cells = group(c0:min(c0+chunkSize-1, numel(group)));
            nChunk = numel(cells);

            % m x nChunk
            localIdx = cellOffset(cells).' + (1:m).';

            dofs = stressDofMap(localIdx);
            sgn = signs(localIdx);

            % --- K_E -> T_E' K_E T_E ---

            % vec(K_E), column-major, m^2 x nChunk
            Kc = reshape(cat(3, K_cons{cells}), m*m, nChunk);
            Ks = reshape(cat(3, K_stab{cells}), m*m, nChunk);

            % vec(s s')
            signOuter = repmat(sgn, m, 1) .* repelem(sgn, m, 1);

            n = m * m * nChunk;
            idx = kPtr + (1:n);

            IK(idx) = repmat(dofs, m, 1);
            JK(idx) = repelem(dofs, m, 1);

            if splitParts
                VKc(idx) = Kc .* signOuter;
                VKs(idx) = Ks .* signOuter;
            else
                VKc(idx) = (Kc + Ks) .* signOuter;
            end

            kPtr = kPtr + n;

            % --- B_E -> B_E T_E ---

            % vec(B_E), column-major, 6m x nChunk
            BE = reshape(cat(3, B_matrix{cells}), 6*m, nChunk);

            % 6(E-1)+i
            dispDofs = 6*(cells(:).' - 1) + (1:6).';

            n = 6 * m * nChunk;
            idx = bPtr + (1:n);

            IB(idx) = repmat(dispDofs, m, 1);
            JB(idx) = repelem(dofs, 6, 1);
            VB(idx) = BE .* repelem(sgn, 6, 1);

            bPtr = bPtr + n;

        end

    end

    %% global sparse assembly

    % sparse() is the peak-memory step
    clear K_cons K_stab B_matrix B_local dofs sgn signOuter Kc Ks BE

    B_global = sparse(IB, JB, VB, nDispDofs, nStressDofs);

    clear IB JB VB

    if splitParts

        K_cons_global = sparse(IK, JK, VKc, nStressDofs, nStressDofs);
        K_stab_global = sparse(IK, JK, VKs, nStressDofs, nStressDofs);

        clear IK JK VKc VKs

        K_global = K_cons_global + K_stab_global;

    else

        K_global = sparse(IK, JK, VKc, nStressDofs, nStressDofs);

    end

end
