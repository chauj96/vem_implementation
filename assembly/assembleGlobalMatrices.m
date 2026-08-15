function [K_cons_global, K_stab_global, B_global] = assembleGlobalMatrices(cell_struct, face_struct, face_global_geom, B_local, K_cons, K_stab)
% assemble global mixed-VEM matrices using a common global face convention
%
% local element contributions are first transformed to the common global
% face convention. matrix entries are collected element-by-element and
% assembled into sparse global matrices only once at the end

    nCells = numel(cell_struct);
    nFaces = numel(face_struct);
    
    nStressDofs = 6 * nFaces;
    nDispDofs = 6 * nCells;
    
    % preallocate assembly storage
    nKEntries = 0;
    nBEntries = 0;
    
    for e = 1:nCells
    
        nLocalStressDofs = 6 * numel(cell_struct(e).faces);
    
        nKEntries = nKEntries + nLocalStressDofs^2;
        nBEntries = nBEntries + 6 * nLocalStressDofs;
    
    end
    
    % K_cons and K_stab have the same sparsity contributions
    IK = zeros(nKEntries,1);
    JK = zeros(nKEntries,1);
    VKc = zeros(nKEntries,1);
    VKs = zeros(nKEntries,1);
    
    % B contribution
    IB = zeros(nBEntries,1);
    JB = zeros(nBEntries,1);
    VB = zeros(nBEntries,1);
    
    kPtr = 0;
    bPtr = 0;
    
    for e = 1:nCells
    
        elementFaces = cell_struct(e).faces;
    
        nElementFaces = numel(elementFaces);
        nLocalStressDofs = 6 * nElementFaces;
    
        % local -> global stress DOF mapping
        stressDofs = zeros(1,nLocalStressDofs);
    
        % element transformation: local face convention -> global face convention
        TE = zeros(nLocalStressDofs);
    
        for lf = 1:nElementFaces
    
            f = elementFaces(lf);
    
            localIdx = 6*(lf-1) + (1:6);
            globalIdx = 6*(f-1) + (1:6);
    
            stressDofs(localIdx) = globalIdx;
    
            localGeom = B_local{e}.geom(lf);
            globalGeom = face_global_geom(f);
    
            Tf = computeFaceTransformation(f, localGeom, globalGeom, face_struct);
    
            assert(norm(Tf' * Tf - eye(6),'fro') < 1e-12, 'Invalid face transformation: cell %d, face %d.', e, f);
    
            TE(localIdx,localIdx) = Tf;
    
        end
    
        % six displacement DOFs belonging to element e
        dispDofs = 6*(e-1) + (1:6);
    
        % local matrices
        Kc = K_cons{e};
        Ks = K_stab{e};
        BE = B_local{e}.matrix;
    
        % transform to global face convention
        Kc = TE' * Kc * TE;
        Ks = TE' * Ks * TE;
        BE = BE * TE;
    
        [rowsK,colsK] = ndgrid(stressDofs,stressDofs);
    
        n = nLocalStressDofs^2;
        idx = kPtr + (1:n);
    
        IK(idx) = rowsK(:);
        JK(idx) = colsK(:);
        VKc(idx) = Kc(:);
        VKs(idx) = Ks(:);
    
        kPtr = kPtr + n;
    
        [rowsB,colsB] = ndgrid(dispDofs,stressDofs);
    
        n = 6 * nLocalStressDofs;
        idx = bPtr + (1:n);
    
        IB(idx) = rowsB(:);
        JB(idx) = colsB(:);
        VB(idx) = BE(:);
    
        bPtr = bPtr + n;
    
    end
    
    % global sparse assembly
    
    K_cons_global = sparse(IK, JK, VKc, nStressDofs, nStressDofs);
    
    K_stab_global = sparse(IK, JK, VKs, nStressDofs, nStressDofs);
    
    B_global = sparse(IB, JB, VB, nDispDofs, nStressDofs);

end
