function [Ared, bred, bc] = applyFlowBC(A, b, cell_struct, face_struct, fracture, flowDofs, exactFun)

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

        assert(numel(localFace) == 1, 'Invalid boundary face connectivity.');

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

    fixedDofs = flowDofs.offset.wF + fracture.edge_flux_dofs(boundaryEdges);

    fixedDofs = unique(fixedDofs(:));

    fixedValues = zeros(numel(fixedDofs), 1);

    %% REDUCED SYSTEM

    freeDofs = setdiff((1:nTotal)', fixedDofs);

    Ared = A(freeDofs, freeDofs);

    bred = b(freeDofs) - A(freeDofs, fixedDofs) * fixedValues;

    %% STORE BC INFORMATION

    bc.freeDofs = freeDofs;
    bc.fixedDofs = fixedDofs;
    bc.fixedValues = fixedValues;
    bc.fullRHS = b;

end
