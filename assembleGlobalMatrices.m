function [K_cons_global, K_stab_global, B_global] = assembleGlobalMatrices(cell_struct, face_struct, face_global_geom, B_local, K_cons, K_stab)
% Assemble global mixed-VEM matrices using a common global face convention.
    
    nCells = numel(cell_struct);
    nFaces = numel(face_struct);
    
    nStressDofs = 6 * nFaces;
    nDispDofs = 6 * nCells;
    
    K_cons_global = sparse(nStressDofs, nStressDofs);
    K_stab_global = sparse(nStressDofs, nStressDofs);
    B_global = sparse(nDispDofs, nStressDofs);
    
    for e = 1:nCells
    
        elementFaces = cell_struct(e).faces;
        nElementFaces = numel(elementFaces);
        nLocalStressDofs = 6 * nElementFaces;
    
        % element transformation: local face convention -> global face convention
        TE = zeros(nLocalStressDofs);
    
        % local -> global stress DOF mapping
        stressDofs = zeros(1, nLocalStressDofs);
    
        for lf = 1:nElementFaces
    
            f = elementFaces(lf);
    
            localIdx = 6*(lf-1) + (1:6);
            globalIdx = 6*(f-1) + (1:6);
    
            stressDofs(localIdx) = globalIdx;
    
            % already-computed local/global face geometries
            localGeom = B_local{e}.geom(lf);
            globalGeom = face_global_geom(f);
    
            % determine local-to-global face transformation
            Tf = computeFaceTransformation(f, localGeom, globalGeom, face_struct);
    
            % basic sanity check
            assert(norm(Tf' * Tf - eye(6), 'fro') < 1e-12, 'Invalid face transformation: cell %d, face %d.', e, f);
    
            TE(localIdx, localIdx) = Tf;
        end
    
        % displacement DOFs belonging to cell e
        dispDofs = 6*(e-1) + (1:6);
    
        % local matrices
        Kc = K_cons{e};
        Ks = K_stab{e};
        BE = B_local{e}.matrix;
    
        % transform to common global face convention
        Kc = TE' * Kc * TE;
        Ks = TE' * Ks * TE;
        BE = BE * TE;
    
        % assemble
        K_cons_global(stressDofs, stressDofs) = K_cons_global(stressDofs, stressDofs) + Kc;
        K_stab_global(stressDofs, stressDofs) = K_stab_global(stressDofs, stressDofs) + Ks;
        B_global(dispDofs, stressDofs) = B_global(dispDofs, stressDofs) + BE;
    
    end

end
