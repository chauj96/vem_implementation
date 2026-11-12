function [Ared, bred, bc] = applyFlowBC( ...
    A, b, cell_struct, face_struct, fracture, flowDofs, exactFun)
% Apply mixed Darcy boundary conditions.
%
% Bulk exterior:
%   p = p_exact (Dirichlet, imposed through flux RHS)
%
% Fracture perimeter:
%   wF . nuF = 0 (Neumann, fixed integrated flux DOFs)
%
% Bulk-fracture Robin coupling is already included in A.

    nTotal = size(A,1);

    %% 1. BULK PRESSURE DIRICHLET BC

    boundaryFaces = find([face_struct.is_pressure_bc]);

    for i = 1:numel(boundaryFaces)

        f = boundaryFaces(i);

        neighbors = face_struct(f).cells(:);
        neighbors = neighbors(neighbors > 0);

        assert(numel(neighbors)==1, ...
            'Pressure BC face must be exterior.');

        e = neighbors(1);

        localFace = find(cell_struct(e).faces == f);

        assert(numel(localFace)==1, ...
            'Invalid boundary face connectivity.');

        s = cell_struct(e).faces_orientation(localFace);

        xf = face_struct(f).center(:);

        exact = exactFun(xf);
        pD = exact.pBulk;

        dof = face_struct(f).flow_dof;

        % Boundary term in the bulk Darcy equation
        b(dof) = b(dof) - s*pD;

    end

    %% 2. FRACTURE PERIMETER NO-FLUX BC

    boundaryEdges = find(fracture.is_flux_bc_edge);

    assert(isequal(sort(boundaryEdges(:)), ...
        sort(fracture.boundary_edges(:))), ...
        'Inconsistent fracture boundary tags.');

    fixedDofs = flowDofs.fractureBoundaryFluxDofs(:);

    fixedValues = ...
        flowDofs.fractureBoundaryFluxValues(:);

    assert(numel(fixedDofs)==numel(fixedValues));

    %% 3. REDUCE SYSTEM

    [fixedDofs, ia] = unique(fixedDofs);
    fixedValues = fixedValues(ia);

    freeDofs = setdiff((1:nTotal)',fixedDofs);

    Ared = A(freeDofs,freeDofs);

    bred = b(freeDofs) ...
         - A(freeDofs,fixedDofs)*fixedValues;

    %% 4. STORE BC INFORMATION

    bc.freeDofs = freeDofs;
    bc.fixedDofs = fixedDofs;
    bc.fixedValues = fixedValues;
    bc.fullRHS = b;

    fprintf('\n[# FLOW BOUNDARY CONDITIONS]\n');
    fprintf('  Bulk pressure BC faces: %d\n', ...
        numel(boundaryFaces));
    fprintf('  Fracture no-flux DOFs:  %d\n', ...
        numel(fixedDofs));

end


% function [Ared, bred, bc] = applyFlowBC(A, b, cell_struct, face_struct, fracture, flowDofs, exactFun)
% 
%     nTotal = size(A, 1);
% 
%     %% BULK EXTERIOR DIRICHLET BC
% 
%     nBoundaryFaces = 0;
% 
%     for f = 1:numel(face_struct)
% 
%         neighbors = face_struct(f).cells(:);
%         neighbors = neighbors(neighbors > 0);
% 
%         if numel(neighbors) ~= 1
%             continue;
%         end
% 
%         e = neighbors(1);
% 
%         localFace = find(cell_struct(e).faces == f);
% 
%         assert(numel(localFace) == 1, 'Invalid boundary face connectivity.');
% 
%         s = cell_struct(e).faces_orientation(localFace);
% 
%         xf = face_struct(f).center(:);
% 
%         exact = exactFun(xf);
%         pD = exact.pBulk;
% 
%         dof = face_struct(f).flow_dof;
% 
%         b(dof) = b(dof) - s * pD;
% 
%         nBoundaryFaces = nBoundaryFaces + 1;
%     end
% 
%     %% FRACTURE BOUNDARY NO-FLUX BC
% 
%     boundaryEdges = fracture.boundary_edges(:);
% 
%     fixedDofs = flowDofs.offset.wF + fracture.edge_flux_dofs(boundaryEdges);
% 
%     fixedDofs = unique(fixedDofs(:));
% 
%     fixedValues = zeros(numel(fixedDofs), 1);
% 
%     %% REDUCED SYSTEM
% 
%     freeDofs = setdiff((1:nTotal)', fixedDofs);
% 
%     Ared = A(freeDofs, freeDofs);
% 
%     bred = b(freeDofs) - A(freeDofs, fixedDofs) * fixedValues;
% 
%     %% STORE BC INFORMATION
% 
%     bc.freeDofs = freeDofs;
%     bc.fixedDofs = fixedDofs;
%     bc.fixedValues = fixedValues;
%     bc.fullRHS = b;
% 
% end
