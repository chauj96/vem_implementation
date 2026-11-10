function b = assembleFlowRHS(cell_struct, fracture, flowDofs, face_struct, exactFun)

    nW  = flowDofs.nW;
    nP  = flowDofs.nP;
    nWF = flowDofs.nWF;
    nPF = flowDofs.nPF;

    bW  = zeros(nW, 1);
    bP  = zeros(nP, 1);
    bWF = zeros(nWF, 1);
    bPF = zeros(nPF, 1);

    %% BULK SOURCE
    assert(numel(cell_struct) == nP, 'Check bulk pressure DOF mapping.');
    
    % two-point Gauss-Legendre quadrature
    g = [-1, 1] / sqrt(3);
    
    for e = 1:numel(cell_struct)
    
        xc = cell_struct(e).center(:);
        Ve = cell_struct(e).volume;
    
        % cartesian cubic mesh
        h = Ve^(1/3);
    
        integralF = 0;
    
        for i = 1:2
            for j = 1:2
                for k = 1:2
    
                    xq = xc + (h/2)*[g(i); g(j); g(k)];
    
                    exact = exactFun(xq);
    
                    integralF = integralF + exact.fBulk;
    
                end
            end
        end
    
        bP(e) = (Ve/8) * integralF;
    end


    %% FRACTURE SOURCE

    for k = 1:fracture.nFaces

        f = fracture.faces(k);
        dof = fracture.pressure_dofs(k);

        xf = face_struct(f).center(:);
        Af = fracture.face_areas(k);

        exact = exactFun(xf);

        bPF(dof) = bPF(dof) + Af * exact.fFrac;
    end

    %% GLOBAL RHS

    b = [bW; bP; bWF; bPF];

end

