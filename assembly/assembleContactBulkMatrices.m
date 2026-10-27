function [A_global, B_global] = assembleContactBulkMatrices(cell_struct, face_struct, face_global_geom, B_local, K_cons, K_stab, nStressDofs)
% assemble bulk mixed VEM operators for elasticity with fracture.
%
% A_global: stress compliance matrix
% B_global: stress-displacement coupling matrix
%
% uses fracture-aware stress DOF numbering stored in cell_struct(e).stress_dofs.

    nCells = numel(cell_struct);
    nDispDofs = 6 * nCells;

    A_global = sparse(nStressDofs, nStressDofs);
    B_global = sparse(nDispDofs, nStressDofs);

    for e = 1:nCells

        elementFaces = cell_struct(e).faces;
        nElementFaces = numel(elementFaces);

        nLocalStressDofs = 6 * nElementFaces;

        % Fracture-aware global stress DOFs
        stressDofs = cell_struct(e).stress_dofs(:).';

        assert(numel(stressDofs) == nLocalStressDofs);

        % Local-to-global transformation
        TE = zeros(nLocalStressDofs);

        for lf = 1:nElementFaces

            f = elementFaces(lf);

            localIdx = 6*(lf-1) + (1:6);

            localGeom = B_local{e}.geom(lf);
            globalGeom = face_global_geom(f);

            Tf = computeFaceTransformation(f, localGeom, globalGeom, face_struct);

            TE(localIdx, localIdx) = Tf;

        end

        % Bulk displacement DOFs
        dispDofs = 6*(e-1) + (1:6);

        % Local mixed VEM matrices
        AE = K_cons{e} + K_stab{e};
        BE = B_local{e}.matrix;

        % Transform to global stress DOF convention
        AE = TE.' * AE * TE;
        BE = BE * TE;

        % Assemble
        A_global(stressDofs, stressDofs) = A_global(stressDofs, stressDofs) + AE;
        B_global(dispDofs, stressDofs) = B_global(dispDofs, stressDofs) + BE;

    end
end
