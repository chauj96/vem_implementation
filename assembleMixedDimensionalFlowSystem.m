function [A, blocks] = assembleMixedDimensionalFlowSystem(cell_struct, fracture, flowDofs, A_mfd_local, E_mfd_local, AF_mfd_local, EF_mfd_local, R, TF)
% assemble the global mixed-dimensional Darcy system.
%
% unknown ordering:
%   x = [w; p; wF; pF]
%
% local MFD inner products use outward-oriented fluxes. global flux DOFs use fixed face/edge orientations.
%
% global matrix:
%
%   [ W+R    E'      0      TF' ]
%   [ -E     0       0       0  ]
%   [  0     0      WF      EF' ]
%   [ -TF    0      -EF      0  ]

    nW = flowDofs.nW;
    nP = flowDofs.nP;
    nWF = flowDofs.nWF;
    nPF = flowDofs.nPF;

    %% 1. Assemble bulk operators

    W = sparse(nW, nW);
    E = sparse(nP, nW);

    for e = 1:numel(cell_struct)

        dofs = cell_struct(e).flow_dofs(:);
        signs = cell_struct(e).flow_signs(:);

        S = diag(signs);

        Wlocal = S * A_mfd_local{e} * S;
        Elocal = E_mfd_local{e};

        pDof = cell_struct(e).pressure_dof;

        W(dofs, dofs) = W(dofs, dofs) + Wlocal;
        E(pDof, dofs) = E(pDof, dofs) + Elocal;
    end

    %% 2. Assemble fracture operators

    WF = sparse(nWF, nWF);
    EF = sparse(nPF, nWF);

    for k = 1:fracture.nFaces

        dofs = fracture.face_flux_dofs{k}(:);
        signs = fracture.face_flux_signs{k}(:);

        S = diag(signs);

        WFlocal = S * AF_mfd_local{k} * S;
        EFlocal = EF_mfd_local{k};

        pFDof = fracture.pressure_dofs(k);

        WF(dofs, dofs) = WF(dofs, dofs) + WFlocal;
        EF(pFDof, dofs) = EF(pFDof, dofs) + EFlocal;
    end

    %% 3. Add interface resistance

    Wstar = W + R;

    %% 4. Assemble global block matrix

    A = [Wstar, E', sparse(nW, nWF), TF';
         -E, sparse(nP, nP), sparse(nP, nWF), sparse(nP, nPF);
         sparse(nWF, nW), sparse(nWF, nP), WF, EF';
         -TF, sparse(nPF, nP), -EF, sparse(nPF, nPF)];

    %% 5. Store assembled blocks

    blocks.W = W;
    blocks.E = E;
    blocks.WF = WF;
    blocks.EF = EF;
    blocks.R = R;
    blocks.TF = TF;
    blocks.Wstar = Wstar;

end
