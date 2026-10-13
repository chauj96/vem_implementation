function [K_global, B_global, D_global, A_global, E_global, M_global] = assembleBiotGlobalMatrices(cell_struct, face_struct, face_global_geom, B_local, K_cons, K_stab, D_biot_local, A_mfd_local, E_mfd_local, M_biot_local)
% assemble global mixed VEM - mixed MFD matrices for the Biot problem
%
% The returned matrices are used in
%
%   [ K    B'    0      D'   ]
%   [-B    0     0      0    ]
%   [ 0    0     A      E'   ]
%   [ D    0   -dt*E    M    ]
%
% No boundary conditions, RHS terms, or time discretization are
% applied in this function.

    nCells = numel(cell_struct);
    nFaces = numel(face_struct);

    nStressDofs = 6 * nFaces;
    nDispDofs = 6 * nCells;

    % initialize global sparse matrices
    K_global = sparse(nStressDofs, nStressDofs);
    B_global = sparse(nDispDofs, nStressDofs);
    D_global = sparse(nCells, nStressDofs);

    A_global = sparse(nFaces, nFaces);
    E_global = sparse(nCells, nFaces);
    M_global = sparse(nCells, nCells);

    for e = 1:nCells

        elementFaces = cell_struct(e).faces;
        nElementFaces = numel(elementFaces);

        %% ========================================================
        % MIXED VEM

        nLocalStressDofs = 6 * nElementFaces;

        % global stress DOFs for this element
        stressDofs = zeros(1, nLocalStressDofs);

        % local -> global stress transformation
        TE = zeros(nLocalStressDofs);

        for lf = 1:nElementFaces

            f = elementFaces(lf);

            localIdx = 6*(lf-1) + (1:6);
            globalIdx = 6*(f-1) + (1:6);

            stressDofs(localIdx) = globalIdx;

            localGeom = B_local{e}.geom(lf);
            globalGeom = face_global_geom(f);

            Tf = computeFaceTransformation(f, localGeom, globalGeom, face_struct);

            TE(localIdx, localIdx) = Tf;

        end

        % displacement DOFs: 6 per cell
        dispDofs = 6*(e-1) + (1:6);

        % local matrices
        KE = K_cons{e} + K_stab{e};
        BE = B_local{e}.matrix;
        DE = D_biot_local{e};

        % transform stress DOFs to common global face convention
        KE = TE' * KE * TE;
        BE = BE * TE;
        DE = DE * TE;

        % assemble VEM matrices
        K_global(stressDofs, stressDofs) = K_global(stressDofs, stressDofs) + KE;

        B_global(dispDofs, stressDofs) = B_global(dispDofs, stressDofs) + BE;

        D_global(e, stressDofs) = D_global(e, stressDofs) + DE;


        %% ========================================================
        % MIXED MFD

        % one global flux DOF per face
        fluxDofs = elementFaces;

        % local MFD inner-product matrix
        AE = A_mfd_local{e};

        % local outward orientation -> global face orientation
        s = cell_struct(e).faces_orientation(:);

        SE = spdiags(s, 0, nElementFaces, nElementFaces);

        AE = SE' * AE * SE;

        % assemble A
        A_global(fluxDofs, fluxDofs) = A_global(fluxDofs, fluxDofs) + AE;

        EE = E_mfd_local{e};

        E_global(e, fluxDofs) = E_global(e, fluxDofs) + EE;

        M_global(e,e) = M_global(e,e) + M_biot_local{e};

    end

end