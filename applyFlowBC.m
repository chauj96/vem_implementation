function [Ared, bred, bc] = applyFlowBC( ...
    A, b, cell_struct, face_struct, fracture, flowDofs, exactFun)

    nTotal = size(A, 1);

    %% BULK EXTERIOR DIRICHLET BC

    nBoundaryFaces = 0;

    for f = 1:numel(face_struct)

        neighbors = face_struct(f).cells(:);
        neighbors = neighbors(neighbors > 0);

        if numel(neighbors) ~= 1
            continue;
        end

        e = neighbors(1);

        localFace = find(cell_struct(e).faces == f);

        assert(numel(localFace) == 1, ...
            'Invalid boundary face connectivity.');

        s = cell_struct(e).faces_orientation(localFace);

        xf = face_struct(f).center(:);

        exact = exactFun(xf);
        pD = exact.pBulk;

        dof = face_struct(f).flow_dof;

        b(dof) = b(dof) - s * pD;

        nBoundaryFaces = nBoundaryFaces + 1;
    end

    %% FRACTURE BOUNDARY NO-FLUX BC

    boundaryEdges = fracture.boundary_edges(:);

    fixedDofs = flowDofs.offset.wF ...
        + fracture.edge_flux_dofs(boundaryEdges);

    fixedDofs = unique(fixedDofs(:));

    fixedValues = zeros(numel(fixedDofs), 1);

    %% REDUCED SYSTEM

    freeDofs = setdiff((1:nTotal)', fixedDofs);

    Ared = A(freeDofs, freeDofs);

    bred = b(freeDofs) ...
        - A(freeDofs, fixedDofs) * fixedValues;

    %% STORE BC INFORMATION

    bc.freeDofs = freeDofs;
    bc.fixedDofs = fixedDofs;
    bc.fixedValues = fixedValues;
    bc.fullRHS = b;

    fprintf('\n[# APPLY FLOW BOUNDARY CONDITIONS]\n');
    fprintf('  Bulk Dirichlet faces: %d\n', nBoundaryFaces);
    fprintf('  Fracture no-flux edges: %d\n', numel(fixedDofs));
    fprintf('  Reduced system: %d x %d\n', size(Ared, 1), size(Ared, 2));
    fprintf('  Status: APPLIED\n');
end


% function [Ared, bred, bc] = applyFlowBC( ...
%     A, b, cell_struct, face_struct, fracture, flowDofs)
% 
%     nTotal = size(A, 1);
% 
%     %% BULK EXTERIOR DIRICHLET BC: p = 1 + abs(x - 0.5)
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
%         assert(numel(localFace) == 1, ...
%             'Invalid boundary face connectivity.');
% 
%         s = cell_struct(e).faces_orientation(localFace);
% 
%         xf = face_struct(f).center(:);
% 
%         pD = 1 + abs(xf(1) - fracture.xFrac);
% 
%         dof = face_struct(f).flow_dof;
% 
%         b(dof) = b(dof) - s * pD;
% 
%         nBoundaryFaces = nBoundaryFaces + 1;
%     end
% 
%     %% FRACTURE BOUNDARY NO-FLUX BC: wF = 0
% 
%     boundaryEdges = fracture.boundary_edges(:);
% 
%     fixedDofs = flowDofs.offset.wF ...
%         + fracture.edge_flux_dofs(boundaryEdges);
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
%     bred = b(freeDofs) ...
%         - A(freeDofs, fixedDofs) * fixedValues;
% 
%     %% STORE BC INFORMATION
% 
%     bc.freeDofs = freeDofs;
%     bc.fixedDofs = fixedDofs;
%     bc.fixedValues = fixedValues;
%     bc.fullRHS = b;
% 
%     fprintf('\n[# APPLY FLOW BOUNDARY CONDITIONS]\n');
%     fprintf('  Bulk Dirichlet faces: %d\n', nBoundaryFaces);
%     fprintf('  Fracture no-flux edges: %d\n', numel(fixedDofs));
%     fprintf('  Reduced system: %d x %d\n', size(Ared, 1), size(Ared, 2));
%     fprintf('  Status: PASS\n');
% end