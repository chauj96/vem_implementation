function fracture = assignInterfaceDofs(fracture)

    nF = fracture.nFaces;

    fracture.traction_dofs = zeros(nF, 6);
    fracture.disp_minus_dofs = zeros(nF, 6);
    fracture.disp_plus_dofs  = zeros(nF, 6);

    for k = 1:nF

        fracture.traction_dofs(k,:) = 6*(k-1) + (1:6);

        fracture.disp_minus_dofs(k,:) = 12*(k-1) + (1:6);

        fracture.disp_plus_dofs(k,:) = 12*(k-1) + (7:12);
    end

    fracture.nTractionDofs = 6*nF;
    fracture.nDispDofs = 12*nF;
end