% clear; clc; close all
addpath('geometry', 'operators', 'assembly', 'solvers', 'io', 'tests', genpath('FACTORIZE'), 'meshes');
% setupMRSTAuto();

%% STEP 1: LOAD MESH
fprintf('\n[# MESH LOADED: INCLINED FRACTURE MODEL]\n');

Lx = 20;
Ly = 20;
Lz = 10;

h = 0.5;
Nz = 1;

psi = 20;
b = 1;
Nfrac = 25;  % number of fracture faces

[G, fractureInfo] = inclinedFracture(Lx, Ly, Lz, h, Nz, psi, b, Nfrac);

[cell_struct, face_struct, V3, cells3D] = MRSTGridConvert(G);

% figure('Color', 'w');
% 
% % plot solid mesh with visible cell edges
% plotGrid(G, 'FaceColor', [0.65, 0.8, 0.95], 'FaceAlpha', 0.35, 'EdgeColor', [0.25, 0.25, 0.25], 'EdgeAlpha', 0.5);
% hold on;
% 
% % highlight fracture faces
% plotFaces(G, fractureInfo.faces, 'FaceColor', [1, 0.2, 0.2], 'FaceAlpha', 0.9, 'EdgeColor', 'k');
% 
% axis equal tight;
% view(3);
% 
% xlabel('x');
% ylabel('y');
% zlabel('z');
% title('3D Inclined Fracture Mesh');
% 
% grid on;
% box on;

%% STEP 2: SET UP GEOMETRY
% Lamé parameters (currently assuming homogeneous isotropic elasticity)
E = 25000;       % MPa
nu = 0.25;

mu = E/(2*(1+nu));
lambda = E*nu/((1+nu)*(1-2*nu));

alpha = 0;       % pure elasticity
s0 = 0;          % no flow/storage in this test

cell_struct = assignCellProperties(cell_struct, V3, cells3D, lambda, mu, alpha, s0);
face_struct = assignFaceProperties(face_struct, Lx, Ly, Lz); 

% expand stress DOFs over fracture
sigma0 = 100;
face_struct = assignContactBoundaryConditions(face_struct, sigma0);
[cell_struct, face_struct, fracture] = assignFractureProperties(cell_struct, face_struct, fractureInfo);

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
cn = 1e6;     % large normal stiffness
cT = 1e-3;    % small tangential stiffness

Ac = computeInterfaceAc(fracture, face_struct, face_global_geom, cn, cT);

J = computeInterfaceJ(fracture, cell_struct, face_struct, B_local, face_global_geom);

%% STEP 4: GLOBAL ASSEMBLY
fprintf('\n[# ASSEMBLE GLOBAL MATRICES]\n');
[A_global, B_global] = assembleContactBulkMatrices(cell_struct, face_struct, face_global_geom, B_local, K_cons, K_stab, nStressDofs);

K_global = assembleContactGlobalSystem(A_global, B_global, Ac, T, J);

RHS = assembleContactRHS(cell_struct, face_struct, face_global_geom, nStressDofs, nDispDofs, nContactTractionDofs, nInterfaceDispDofs);
[K_bc, RHS_bc, bc] = applyContactBoundaryConstraints(K_global, RHS, cell_struct, face_struct, face_global_geom);

%% STEP 5: SOLVE GLOBAL CONTACT SYSTEM
fprintf('\n[# SOLVE CONTACT SYSTEM]\n');

y = K_bc \ RHS_bc;
solution = y(1:bc.nOriginal);

%% STEP 6: EXTRACT SOLUTION
fprintf('\n[# EXTRACT SOLUTION]\n');

i1 = nStressDofs;
i2 = i1 + nDispDofs;
i3 = i2 + nContactTractionDofs;

sigma_h = solution(1:i1);
u_h = solution(i1+1:i2);
tc_h = solution(i2+1:i3);
uc_h = solution(i3+1:end);

%% STEP 7: PLOT FRACTURE TRACTION
nF = fracture.nFaces;

jump_N = zeros(nF,1);
jump_T = zeros(nF,1);
xF = zeros(nF,1);
yF = zeros(nF,1);

for k = 1:nF

    f = fracture.faces(k);

    xf = face_struct(f).center(:);
    Af = face_struct(f).area;

    gc = face_global_geom(f);

    % minus and plus cells
    em = fracture.minus_cells(k);
    ep = fracture.plus_cells(k);

    % local face indices
    lfm = find(cell_struct(em).faces == f);
    lfp = find(cell_struct(ep).faces == f);

    % local face geometries
    gm = B_local{em}.geom(lfm);
    gp = B_local{ep}.geom(lfp);

    % mortar displacement DOFs
    um = uc_h(fracture.disp_minus_dofs(k,:));
    up = uc_h(fracture.disp_plus_dofs(k,:));

    % evaluate local displacement basis at face center
    Phi_m = evaluateFaceBasis(xf, xf, gm.Qf, gm.t1, gm.t2, gm.n, Af, gm.m20, gm.m02);

    Phi_p = evaluateFaceBasis(xf, xf, gp.Qf, gp.t1, gp.t2, gp.n, Af, gp.m20, gp.m02);

    % reconstruct physical displacements on each side
    u_minus = Phi_m * um;
    u_plus  = Phi_p * up;
    
    % displacement jump: plus minus minus
    jump = u_plus - u_minus;
    
    % project onto the fixed fracture directions
    jump_N(k) = dot(jump, fracture.normal);
    jump_T(k) = dot(jump, fracture.tangent);

    xF(k) = xf(1);
    yF(k) = xf(2);

end

%% STEP 7: NUMERICAL VS ANALYTICAL FRACTURE JUMPS

% position along fracture
s = xF*cosd(psi) + yF*sind(psi);
[s, idx] = sort(s);

jump_N_num = jump_N(idx);
jump_T_num = jump_T(idx);

% analytical solution: frictionless closed fracture (plane strain)

s_exact = linspace(-b, b, 500);
tau0 = sigma0 * sind(psi) * cosd(psi);
jump_N_exact = zeros(size(s_exact));

% analytical tangential jump
C = 4 * (1-nu^2) * tau0/E;

jump_T_exact = C * sqrt(max(b^2 - s_exact.^2, 0));

% evaluate analytical solutions at numerical face centers
jump_N_exact_faces = zeros(size(s));
jump_T_exact_faces = C * sqrt(max(b^2 - s.^2, 0));


figure('Color','w','Position',[200 150 850 650]);


subplot(2,1,1)

plot(s, jump_N_num, 'bo-', 'LineWidth',1.5,'MarkerSize',5);
hold on;

plot(s_exact, jump_N_exact, 'r--', 'LineWidth',2);

xlabel('Position along fracture, s (m)');
ylabel('Normal displacement jump (m)');
title('Normal Displacement Jump');
legend('Mixed VEM','Analytical','Location','best');
grid on;
xlim([-b b]);

subplot(2,1,2)

plot(s, jump_T_num, 'bo-', 'LineWidth',1.5,'MarkerSize',5);
hold on;

plot(s_exact, jump_T_exact, 'r--', 'LineWidth',2);

xlabel('Position along fracture, s (m)');
ylabel('Tangential displacement jump (m)');
title('Tangential Displacement Jump');
legend('Mixed VEM','Analytical','Location','best');
grid on;
xlim([-b b]);

sgtitle('Inclined Fracture under Compression');


%% STEP 8: CONTACT TRACTION RECONSTRUCTION

nF = fracture.nFaces;

traction_N = zeros(nF,1);
traction_T = zeros(nF,1);

nGamma = fracture.normal(:);
tGamma = fracture.tangent(:);

for k = 1:nF

    f = fracture.faces(k);

    xf = face_struct(f).center(:);
    Af = face_struct(f).area;

    gc = face_global_geom(f);

    % contact traction DOFs
    dofs = fracture.traction_dofs(k,:);
    tc_local = tc_h(dofs);

    % reconstruct contact traction at face center
    Phi_c = evaluateFaceBasis(xf, xf, gc.Qf, gc.t1, gc.t2, gc.n, Af, gc.m20, gc.m02);

    tc = Phi_c * tc_local;

    % project onto fixed fracture directions
    traction_N(k) = dot(tc, nGamma);
    traction_T(k) = dot(tc, tGamma);

end

traction_N_num = traction_N(idx);
traction_T_num = traction_T(idx);

% analytical traction

traction_N_exact = -sigma0*sind(psi)^2 * ones(size(s_exact));
traction_T_exact = zeros(size(s_exact));

traction_T_reference = cT * jump_T_exact;
traction_T_from_jump = cT * jump_T_num;

figure('Color','w','Position',[250 180 850 650]);

% normal traction
subplot(2,1,1);

plot(s, traction_N_num, 'bo-', 'LineWidth',1.5,'MarkerSize',5);
hold on;

plot(s_exact, traction_N_exact, 'r--', 'LineWidth',2);

xlabel('Position along fracture, s (m)');
ylabel('Normal traction (MPa)');
title('Normal Contact Traction');

legend('Mixed VEM','Analytical', 'Location','best');

grid on;
xlim([-b b]);

% tangential traction
subplot(2,1,2);

plot(s, traction_T_num, 'bo-', 'LineWidth',1.5,'MarkerSize',5);
hold on;

plot(s_exact, traction_T_exact, 'k-', 'LineWidth',1.5);
plot(s_exact, traction_T_reference, 'r--', 'LineWidth',2);
plot(s, traction_T_from_jump, 'ms:', 'LineWidth',1.2,'MarkerSize',4);

xlabel('Position along fracture, s (m)');
ylabel('Tangential traction (MPa)');
title('Tangential Contact Traction');

legend('Mixed VEM', 'Exact frictionless', 'c_T g_T^{exact}', 'c_T g_{T,h}', 'Location','best');

grid on;
xlim([-b b]);

sgtitle('Inclined Fracture under Compression');
