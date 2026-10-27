function [cell_struct, face_struct, nStressDofs] = assignFractureStressDofs(cell_struct, face_struct, fracture)

    nFaces = numel(face_struct);
    nCells = numel(cell_struct);

    % base stress DOFs: 6 per geometric face
    for f = 1:nFaces
        face_struct(f).stress_dofs = (6*(f-1)+1):(6*f);
        face_struct(f).stress_dofs_plus = [];
    end

    nextDof = 6*nFaces + 1;

    % duplicate stress DOFs on fracture faces
    for k = 1:fracture.nFaces

        f = fracture.faces(k);

        face_struct(f).stress_dofs_plus = nextDof:(nextDof+5);

        nextDof = nextDof + 6;
    end

    % build cell-local to global stress DOF mapping
    for c = 1:nCells

        faces = cell_struct(c).faces;
        localDofs = zeros(1, 6*numel(faces));

        for lf = 1:numel(faces)

            f = faces(lf);

            if face_struct(f).is_fracture && c == fracture.plus_cells(face_struct(f).fracture_id)

                dofs = face_struct(f).stress_dofs_plus;

            else

                dofs = face_struct(f).stress_dofs;

            end

            localDofs(6*(lf-1)+(1:6)) = dofs;
        end

        cell_struct(c).stress_dofs = localDofs;
    end

    nStressDofs = nextDof - 1;

    fprintf('\n[# STRESS DOF NUMBERING]\n');
    fprintf('  Total stress DOFs: %d\n', nStressDofs);
end