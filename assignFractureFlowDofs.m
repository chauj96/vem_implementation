function [cell_struct, face_struct, fracture, flowDofs] = assignFractureFlowDofs(cell_struct, face_struct, fracture)
% assign DOFs for mixed-dimensional Darcy flow.
%
% unknowns:
%   w  : bulk face normal flux
%   p  : bulk cell pressure
%   wF : fracture edge co-normal flux
%   pF : fracture face pressure
%
% each unknown uses block-local indexing. global offsets are stored in flowDofs.offset.

    nFaces = numel(face_struct);
    nCells = numel(cell_struct);

    %% 1. Bulk flux DOFs: w

    % one DOF per geometric face
    for f = 1:nFaces
        face_struct(f).flow_dof = f;
        face_struct(f).flow_dof_plus = [];
    end

    nextDof = nFaces + 1;

    % duplicate flux DOFs on fracture faces
    for k = 1:fracture.nFaces

        f = fracture.faces(k);

        face_struct(f).flow_dof_plus = nextDof;
        nextDof = nextDof + 1;
    end

    flowDofs.nW = nextDof - 1;

    % cell-local flux DOF mapping
    for c = 1:nCells

        faces = cell_struct(c).faces;
        nLocalFaces = numel(faces);

        localDofs = zeros(1, nLocalFaces);

        for lf = 1:nLocalFaces

            f = faces(lf);

            if face_struct(f).is_fracture && c == fracture.plus_cells(face_struct(f).fracture_id)

                localDofs(lf) = face_struct(f).flow_dof_plus;

            else

                localDofs(lf) = face_struct(f).flow_dof;
            end
        end

        cell_struct(c).flow_dofs = localDofs;

        % global normal -> local outward normal
        cell_struct(c).flow_signs = cell_struct(c).faces_orientation;
    end

    %% 2. Bulk pressure DOFs: p

    for c = 1:nCells
        cell_struct(c).pressure_dof = c;
    end

    flowDofs.nP = nCells;

    %% 3. Fracture tangential flux DOFs: wF

    fracture.edge_flux_dofs = (1:fracture.nEdges)';

    fracture.face_flux_dofs = cell(fracture.nFaces, 1);
    fracture.face_flux_signs = fracture.face_edge_signs;

    for k = 1:fracture.nFaces

        fracture.face_flux_dofs{k} = fracture.edge_flux_dofs(fracture.face_edges{k})';
    end

    flowDofs.nWF = fracture.nEdges;

    %% 4. Fracture pressure DOFs: pF

    fracture.pressure_dofs = (1:fracture.nFaces)';

    flowDofs.nPF = fracture.nFaces;

    %% 5. Global block offsets

    flowDofs.offset.w = 0;

    flowDofs.offset.p = flowDofs.nW;

    flowDofs.offset.wF = flowDofs.nW + flowDofs.nP;

    flowDofs.offset.pF = ...
        flowDofs.nW + flowDofs.nP + flowDofs.nWF;

    flowDofs.total = flowDofs.nW + flowDofs.nP + flowDofs.nWF + flowDofs.nPF;

    %% 6. Fracture boundary DOF information

    fracture.boundary_flux_dofs = fracture.edge_flux_dofs(fracture.boundary_edges);

    fracture.internal_flux_dofs = fracture.edge_flux_dofs(fracture.internal_edges);

    fprintf('\n[# FLOW DOF NUMBERING]\n');
    fprintf('  Bulk flux w:       %d\n', flowDofs.nW);
    fprintf('  Bulk pressure p:   %d\n', flowDofs.nP);
    fprintf('  Fracture flux wF:  %d\n', flowDofs.nWF);
    fprintf('  Fracture pressure: %d\n', flowDofs.nPF);
    fprintf('  Total DOFs:        %d\n', flowDofs.total);
end
