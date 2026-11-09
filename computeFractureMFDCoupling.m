function EF_mfd_local = computeFractureMFDCoupling(fracture)
% construct local fracture negative-divergence operators.
%
% E_F,k = -[s_1, ..., s_ne]
%
% signs map globally oriented fracture edge fluxes
% to local outward co-normal fluxes.

    nFaces = fracture.nFaces;

    EF_mfd_local = cell(nFaces, 1);

    for k = 1:nFaces

        s = fracture.face_edge_signs{k}(:)';

        EF_mfd_local{k} = -s;
    end
end
