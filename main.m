% clear; clc; close all
addpath('geometry', 'operators', 'assembly', 'solvers', 'io', 'tests', genpath('FACTORIZE'), 'meshes');
% setupMRSTAuto();

%% STEP 1: LOAD MESH AND CHOOSE SYSTEM (mixed or hybrid), SOLVER TYPE (direct or iterative)
[meshOption, solverOption] = deal(1, 'iterative');

if meshOption == 1
    % unit cube

    fprintf('\n[# MESH LOADED: unit cube model]\n');

    meshName = 'unit_cube';

    N = 32;
    G = cartGrid([N N N],[1 1 1]);
    G = computeGeometry(G);
    [cell_struct, face_struct, V3, cells3D] = MRSTGridConvert(G);

elseif meshOption == 2
    % two-fault model (need to use iterative solver)

    fprintf('\n[# MESH LOADED: two-fault model]\n');

    meshName = 'two_fault';

    % option 1: read .vtu file
    filename = 'meshes/two_fault/two_fault.vtu';
    [cell_struct, face_struct, V3, cells3D] = readVTU(filename);

    % option 2: load preprocessed mesh file
    % load('meshes/two_fault/two_fault.mat'); 

elseif meshOption == 3
    % polyhedral voronoi model (either direct or iterative solver is fine)

    fprintf('\n[# MESH LOADED: polyhedral voronoi model]\n');

    meshName = 'polyhedral_voronoi';

    % option 1: read .vtu file
    filename = 'meshes/polyhedral_voronoi/polyhedral_voronoi_complex.vtu';
    [cell_struct, face_struct, V3, cells3D] = readVTU(filename);

    % option 2: load proprocessed mesh file
    % load('meshes/polyhedral_voronoi/polyhedral_voronoi_complex.mat');

end

%% STEP 2: SET UP GEOMETRY
% Lamé parameters (currently assuming homogeneous isotropic elasticity)
lambda = 1;
mu = 1;
alpha = 1;
s0 = 1;

cell_struct = assignCellProperties(cell_struct, V3, cells3D, lambda, mu, alpha, s0);
face_struct = assignFaceProperties(face_struct); 

cell_struct = computeVolumeQuadrature(cell_struct, face_struct, V3);
face_struct = computeFaceQuadrature(face_struct, V3);

nStressDofs = 6*numel(face_struct);
nDispDofs = 6*numel(cell_struct);

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

% BIOT STRESS-PRESSURE COUPLING
D_biot_local = computeBiotCoupling(cell_struct, P_local);

K_cons = computeConsistency(cell_struct, P_local);
K_stab = computeStabilization(cell_struct, face_struct, B_local, P_local);

% FLOW PART (MIXED MFD)
A_mfd_local = computeMFDInnerProduct(cell_struct, face_struct);
E_mfd_local = computeMFDCoupling(cell_struct);
M_biot_local = computeBiotStorage(cell_struct);

%% STEP 4: GLOBAL ASSEMBLY
fprintf('\n[# ASSEMBLE GLOBAL MATRICES]\n');

[K_global, B_global, D_global, A_global, E_global, M_global] = assembleBiotGlobalMatrices( ...
        cell_struct, face_struct, face_global_geom, ...
        B_local, K_cons, K_stab, D_biot_local, ...
        A_mfd_local, E_mfd_local, M_biot_local);

%% STEP 5: BUILD FULL BIOT SYSTEM
dt = 0.5;

nStressDofs = size(K_global, 1);
nDispDofs = size(B_global, 1);
nFluxDofs = size(A_global, 1);
nPressureDofs = size(M_global, 1);

LHS = [ ...
    K_global,   B_global',   sparse(nStressDofs, nFluxDofs),  D_global';
   -B_global,   sparse(nDispDofs, nDispDofs), sparse(nDispDofs, nFluxDofs), sparse(nDispDofs, nPressureDofs);
    sparse(nFluxDofs, nStressDofs), sparse(nFluxDofs, nDispDofs), A_global, E_global';
    D_global,   sparse(nPressureDofs, nDispDofs), -dt * E_global, M_global ];

[f_sigma, f_u, f_w, f_p] = assembleBiotRHS(cell_struct, face_struct, face_global_geom);

[sigma_old, p_old] = initializeBiot(cell_struct, face_struct, face_global_geom);

RHS = [ ...
       f_sigma;
       f_u;
       f_w;
       D_global * sigma_old + M_global * p_old + dt * f_p ];

solution = LHS \ RHS;

%% STEP 6: SOLUTION CHECK
i1 = nStressDofs;
i2 = i1 + nDispDofs;
i3 = i2 + nFluxDofs;

sigma_h = solution(1:i1);
u_h = solution(i1+1:i2);
w_h = solution(i2+1:i3);
p_h = solution(i3+1:end);

problem.u = @(x) [ ...
    x(1,:).^2;
    x(2,:).^2;
    x(3,:).^2 ];

problem.p = @(x) ...
    x(1,:).^2 + x(2,:).^2 + x(3,:).^2;

problem.w = @(x) [ ...
    -2*x(1,:);
    -2*x(2,:);
    -2*x(3,:) ];

problem.divsigma = @(x) [ ...
    6 - 2*x(1,:);
    6 - 2*x(2,:);
    6 - 2*x(3,:) ];

problem.sigma = @exactSigma;


err = computeBiotErrors(cell_struct, face_struct, face_global_geom, V3, B_local, D_local, P_local, A_mfd_local, sigma_h, u_h, w_h, p_h, problem);

fprintf('\n=====================================================\n');
fprintf('BIOT SPATIAL ERRORS\n');
fprintf('=====================================================\n');
fprintf('Displacement error       : %.6e\n', err.displacement);
fprintf('Stress divergence error  : %.6e\n', err.divergence);
fprintf('Stress projection error  : %.6e\n', err.projection);
fprintf('Stress traction error    : %.6e\n', err.stress);
fprintf('Pressure error           : %.6e\n', err.pressure);
fprintf('Flux error               : %.6e\n', err.flux);
fprintf('=====================================================\n');




function S = exactSigma(x)

    nq = size(x,2);
    S = zeros(3,3,nq);

    for q = 1:nq

        X = x(1,q);
        Y = x(2,q);
        Z = x(3,q);

        p = X^2 + Y^2 + Z^2;

        S(:,:,q) = [ ...
            6*X + 2*Y + 2*Z - p, 0,                       0;
            0,                       2*X + 6*Y + 2*Z - p, 0;
            0,                       0,                       2*X + 2*Y + 6*Z - p ];

    end

end


% fprintf('\n=====================================================\n');
% fprintf('FULL BIOT SOLUTION CHECK\n');
% fprintf('=====================================================\n');
% 
% % Extract solution blocks
% i1 = nStressDofs;
% i2 = i1 + nDispDofs;
% i3 = i2 + nFluxDofs;
% 
% sigma_h = solution(1:i1);
% u_h     = solution(i1+1:i2);
% w_h     = solution(i2+1:i3);
% p_h     = solution(i3+1:end);
% 
% % Print solution
% fprintf('\nsigma_h =\n');
% disp(sigma_h);
% 
% fprintf('u_h =\n');
% disp(u_h);
% 
% fprintf('w_h =\n');
% disp(w_h);
% 
% fprintf('p_h =\n');
% disp(p_h);
% 
% % Residual check
% residual = LHS * solution - RHS;
% 
% fprintf('-----------------------------------------------------\n');
% fprintf('Linear system residual\n');
% fprintf('||LHS*x - RHS||_2       = %.3e\n', norm(residual));
% fprintf('Relative residual       = %.3e\n', ...
%     norm(residual) / max(norm(RHS), eps));
% 
% fprintf('=====================================================\n');

