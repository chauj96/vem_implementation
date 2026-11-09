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
b = 6;

[G, fractureInfo] = inclinedFracture(Lx, Ly, Lz, h, Nz, psi, b);
[cell_struct, face_struct, V3, cells3D] = MRSTGridConvert(G);

figure('Color', 'w');

% plot solid mesh with visible cell edges
plotGrid(G, 'FaceColor', [0.65, 0.8, 0.95], 'FaceAlpha', 0.35, 'EdgeColor', [0.25, 0.25, 0.25], 'EdgeAlpha', 0.5);
hold on;

% highlight fracture faces
plotFaces(G, fractureInfo.faces, 'FaceColor', [1, 0.2, 0.2], 'FaceAlpha', 0.9, 'EdgeColor', 'k');

axis equal tight;
view(3);

xlabel('x');
ylabel('y');
zlabel('z');
title('3D Inclined Fracture Mesh');

grid on;
box on;

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
solution = bc.xp + bc.Z*y;

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

xi = zeros(nF,1);
tN = zeros(nF,1);
tT = zeros(nF,1);

n = fracture.normal(:);
t = fracture.tangent(:);

centers = zeros(3,nF);

for k = 1:nF
    centers(:,k) = face_struct(fracture.faces(k)).center(:);
end

% Coordinate measured from the first fracture face
s = t.' * centers;
s0 = min(s);

for k = 1:nF

    f = fracture.faces(k);
    xf = face_struct(f).center(:);
    geom = face_global_geom(f);

    xi(k) = dot(xf,t) - s0;

    dofs = fracture.traction_dofs(k,:);
    coeff = tc_h(dofs);

    Psi = evaluateFaceBasis(xf, xf, geom.Qf, geom.t1, geom.t2, geom.n, face_struct(f).area, geom.m20, geom.m02);

    traction = Psi * coeff;

    tN(k) = dot(traction,n);
    tT(k) = dot(traction,t);

end

[xi,idx] = sort(xi);
tN = tN(idx);
tT = tT(idx);

figure;

subplot(2,1,1);
plot(xi,tN,'o-','LineWidth',1.5);
yline(-sigma0*sind(psi)^2,'--');
xlabel('\xi');
ylabel('Normal traction (MPa)');
title('Normal Contact Traction');
grid on;

subplot(2,1,2);
plot(xi,tT,'o-','LineWidth',1.5);
xlabel('\xi');
ylabel('Tangential traction (MPa)');
title('Tangential Contact Traction');
grid on;