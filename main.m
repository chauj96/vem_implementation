% clear; clc; close all
addpath('geometry', 'operators', 'assembly', 'solvers', 'io', 'tests', genpath('FACTORIZE'), 'meshes');
% setupMRSTAuto();

%% STEP 1: LOAD MESH
fprintf('\n[# MESH LOADED: unit cube model]\n');

meshName = 'unit_cube';

Nx = 8;
Ny = 8;
Nz = 4;
G = cartGrid([Nx Ny Nz],[1 1 1]);
G = computeGeometry(G);
[cell_struct, face_struct, V3, cells3D] = MRSTGridConvert(G);

%% STEP 2: SET UP GEOMETRY
% Lamé parameters (currently assuming homogeneous isotropic elasticity)
lambda = 1;
mu = 1;
alpha = 1;
s0 = 1;
xFrac = 0.5;

cell_struct = assignCellProperties(cell_struct, V3, cells3D, lambda, mu, alpha, s0);
face_struct = assignFaceProperties(face_struct); 

% expand stress DOFs over fracture (sigma_h)
[cell_struct, face_struct, fracture] = assignFractureProperties(cell_struct, face_struct, xFrac);

cell_struct = computeVolumeQuadrature(cell_struct, face_struct, V3);
face_struct = computeFaceQuadrature(face_struct, V3);

[cell_struct, face_struct, nStressDofs] = assignFractureStressDofs(cell_struct, face_struct, fracture);

% bulk displacement DOFs
nDispDofs = 6*numel(cell_struct);

fracture = assignInterfaceDofs(fracture);

nContactTractionDofs = fracture.nTractionDofs;
nInterfaceDispDofs = fracture.nDispDofs;

nTotalDofs = nStressDofs + nDispDofs + nContactTractionDofs + nInterfaceDispDofs;

%% STEP 3: (LOCAL LEVEL) COMPUTE VEM MATRICES
fprintf('\n[# COMPUTE LOCAL OPERATORS]\n');

% INERTIA MATRIX (quadrature approach, Tonon approach)
A_local = computeAE_Tonon(cell_struct, face_struct, V3);

% DIVERGENCE MATRIX (quadrature approach, closed form approach from notes)
B_local = computeBE(cell_struct, face_struct, V3);
face_global_geom = computeGlobalFaceGeometry(face_struct, V3);

% PROJECTION MATRIX
D_local = computeDE(cell_struct, A_local, B_local);
Bproj = computeBproj(cell_struct, face_struct, D_local, B_local);
P_local = computeProjectionMatrix(cell_struct, Bproj);

K_cons = computeConsistency(cell_struct, P_local);
K_stab = computeStabilization(cell_struct, face_struct, B_local, P_local);

T = computeInterfaceT(fracture, cell_struct, face_struct, B_local, face_global_geom, nStressDofs);

% contact stiffness parameters
cn = 1;
cT = 1;

Ac = computeInterfaceAc(fracture, face_struct, face_global_geom, cn, cT);

J = computeInterfaceJ(fracture, cell_struct, face_struct, B_local, face_global_geom);

%% STEP 4: GLOBAL ASSEMBLY
fprintf('\n[# ASSEMBLE GLOBAL MATRICES]\n');
[A_global, B_global] = assembleContactBulkMatrices(cell_struct, face_struct, face_global_geom, B_local, K_cons, K_stab, nStressDofs);

K_global = assembleContactGlobalSystem(A_global, B_global, Ac, T, J);

RHS = assembleContactRHS(cell_struct, face_struct, face_global_geom, nStressDofs, nDispDofs, nContactTractionDofs, nInterfaceDispDofs, xFrac);

%% STEP 5: SOLVE GLOBAL CONTACT SYSTEM
fprintf('\n[# SOLVE CONTACT SYSTEM]\n');

solution = K_global \ RHS;

% extract unknowns
i1 = nStressDofs;
i2 = i1 + nDispDofs;
i3 = i2 + nContactTractionDofs;

sigma_h = solution(1:i1);
u_h = solution(i1+1:i2);
tc_h = solution(i2+1:i3);
uc_h = solution(i3+1:end);

results = checkContactPatchTest(cell_struct, face_struct, fracture, face_global_geom, B_local, sigma_h, u_h, tc_h, uc_h, A_global, B_global, Ac, T, J, RHS, xFrac);


function results = checkContactPatchTest(cell_struct, face_struct, fracture, face_global_geom, B_local, sigma_h, u_h, tc_h, uc_h, A_global, B_global, Ac, T, J, RHS, xFrac)
    fprintf('\n========================================\n');
    fprintf('       CONTACT PATCH TEST CHECK\n');
    fprintf('========================================\n');

    % sigmaExact = diag([1,0,0]);
    sigmaExact = diag([1, 1/3, 1/3]);
    tExact = [1;0;0];

    %% 1. BULK STRESS
    sigma_exact_dofs = zeros(size(sigma_h));
    visited = false(size(sigma_h));

    for e = 1:numel(cell_struct)

        faces = cell_struct(e).faces;

        for lf = 1:numel(faces)

            f = faces(lf);
            geom = face_global_geom(f);
            Af = face_struct(f).area;

            traction = sigmaExact * geom.n;

            coeff = zeros(6,1);
            coeff(1:3) = Af * [ ...
                dot(traction,geom.t1);
                dot(traction,geom.t2);
                dot(traction,geom.n)];

            cols = 6*(lf-1)+(1:6);
            dofs = cell_struct(e).stress_dofs(cols);

            if any(visited(dofs))
                assert(norm(sigma_exact_dofs(dofs)-coeff) < 1e-10, 'Inconsistent shared stress DOF values.');
            end

            sigma_exact_dofs(dofs) = coeff;
            visited(dofs) = true;

        end
    end

    assert(all(visited), 'Some stress DOFs were not visited.');

    results.stress_error = norm(sigma_h-sigma_exact_dofs) / max(norm(sigma_exact_dofs),eps);

    fprintf('\n[1] BULK STRESS\n');
    fprintf('  Relative DOF error: %.6e\n', results.stress_error);


    %% 2. BULK DISPLACEMENT
    % displacement evaluated at the cell barycenter.

    u_exact_centers = zeros(3*numel(cell_struct),1);
    u_num_centers = zeros(3*numel(cell_struct),1);

    for e = 1:numel(cell_struct)

        xE = cell_struct(e).center(:);
        isPlus = xE(1) > xFrac;
        ue = [xE(1)/3 + double(isPlus); 0; 0];
 
        % 
        % ue = [0.4*xE(1);
        %      -0.1*xE(2);
        %      -0.1*xE(3)];

        % if isPlus
        %     ue(1) = ue(1)+1;
        % end

        u_exact_centers(3*(e-1)+(1:3)) = ue;
        u_num_centers(3*(e-1)+(1:3)) = u_h(6*(e-1)+(1:3));

    end

    results.displacement_center_error = norm(u_num_centers-u_exact_centers) / max(norm(u_exact_centers),eps);

    fprintf('\n[2] BULK DISPLACEMENT\n');
    fprintf('  Relative cell-center error: %.6e\n', results.displacement_center_error);


    %% 3. CONTACT TRACTION
    tc_exact_dofs = zeros(size(tc_h));

    for k = 1:fracture.nFaces

        f = fracture.faces(k);
        geom = face_global_geom(f);
        Af = face_struct(f).area;

        coeff = zeros(6,1);
        coeff(1:3) = Af * [ ...
            dot(tExact,geom.t1);
            dot(tExact,geom.t2);
            dot(tExact,geom.n)];

        dofs = fracture.traction_dofs(k,:);
        tc_exact_dofs(dofs) = coeff;

    end

    results.contact_traction_error = norm(tc_h-tc_exact_dofs) / max(norm(tc_exact_dofs),eps);

    fprintf('\n[3] CONTACT TRACTION\n');
    fprintf('  Relative DOF error: %.6e\n', results.contact_traction_error);


    %% 4. INTERFACE DISPLACEMENT
    % reconstruct displacement on each side using the corresponding local face basis.

    errMinus2 = 0;
    errPlus2 = 0;
    errJump2 = 0;

    normMinus2 = 0;
    normPlus2 = 0;
    normJump2 = 0;

    for k = 1:fracture.nFaces

        f = fracture.faces(k);
        face = face_struct(f);

        em = fracture.minus_cells(k);
        ep = fracture.plus_cells(k);

        lfm = find(cell_struct(em).faces == f);
        lfp = find(cell_struct(ep).faces == f);

        gm = B_local{em}.geom(lfm);
        gp = B_local{ep}.geom(lfp);

        dm = fracture.disp_minus_dofs(k,:);
        dp = fracture.disp_plus_dofs(k,:);

        cm = uc_h(dm);
        cp = uc_h(dp);

        xf = face.center(:);
        Af = face.area;

        qp = face.quad_points;
        w = face.quad_weights;

        for q = 1:numel(w)

            xq = qp(:,q);

            PhiMinus = evaluateFaceBasis(xq,xf,gm.Qf,gm.t1,gm.t2,gm.n, Af,gm.m20,gm.m02);

            PhiPlus = evaluateFaceBasis(xq,xf,gp.Qf,gp.t1,gp.t2,gp.n, Af,gp.m20,gp.m02);

            um_num = PhiMinus*cm;
            up_num = PhiPlus*cp;

            % um_exact = [0.4*xq(1);
            %            -0.1*xq(2);
            %            -0.1*xq(3)];
            % 
            % up_exact = um_exact+[1;0;0];

            um_exact = [xq(1)/3; 0; 0];
            up_exact = um_exact + [1; 0; 0];

            jump_num = up_num-um_num;
            jump_exact = up_exact-um_exact;

            errMinus2 = errMinus2 + w(q)*sum((um_num-um_exact).^2);

            errPlus2 = errPlus2 + w(q)*sum((up_num-up_exact).^2);

            errJump2 = errJump2 + w(q)*sum((jump_num-jump_exact).^2);

            normMinus2 = normMinus2 + w(q)*sum(um_exact.^2);

            normPlus2 = normPlus2 + w(q)*sum(up_exact.^2);

            normJump2 = normJump2 + w(q)*sum(jump_exact.^2);

        end
    end

    results.interface_minus_error = sqrt(errMinus2/max(normMinus2,eps));

    results.interface_plus_error = sqrt(errPlus2/max(normPlus2,eps));

    results.jump_error = sqrt(errJump2/max(normJump2,eps));

    fprintf('\n[4] INTERFACE DISPLACEMENT\n');
    fprintf('  Minus-side relative error: %.6e\n', results.interface_minus_error);
    fprintf('  Plus-side relative error:  %.6e\n', results.interface_plus_error);
    fprintf('  Jump relative error:       %.6e\n', results.jump_error);

end
