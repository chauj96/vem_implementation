% clear; clc; close all
addpath('geometry', 'operators', 'assembly', 'solvers', 'io', 'tests', genpath('FACTORIZE'));
% setupMRSTAuto();  % for MRST automatic setup

%% STEP 1: LOAD MESH AND CHOOSE SOLVER TYPE (direct or iterative)
meshOption = 3;
solverOption = 'iterative'; 

if meshOption == 1
    % unit cube

    fprintf('\n[# MESH LOADED: unit cube model]\n');

    meshName = 'unit_cube';

    G = cartGrid([10 10 10],[1 1 1]);
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

cell_struct = assignCellProperties(cell_struct, V3, cells3D, lambda, mu);
face_struct = assignFaceProperties(face_struct); 

cell_struct = computeVolumeQuadrature(cell_struct, face_struct, V3);
face_struct = computeFaceQuadrature(face_struct, V3);

nStressDofs = 6*numel(face_struct);
nDispDofs = 6*numel(cell_struct);

%% STEP 3: COMPUTE VEM MATRICES (LOCAL LEVEL)
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

%% STEP 4: ASSEMBLE GLOBAL SYSTEM (GLOBAL LEVEL)
fprintf('\n[# ASSEMBLE GLOBAL SYSTEM]\n');

[K_cons_global, K_stab_global, B_global] = assembleGlobalMatrices(cell_struct, face_struct, face_global_geom, B_local, K_cons, K_stab);

K_global = K_cons_global + K_stab_global;
A_global = [K_global, B_global';
            B_global, sparse(6*numel(cell_struct), 6*numel(cell_struct))];

rhs_D = assembleDirichletRHS(cell_struct, face_struct, face_global_geom);

% f = 0 (just for now)
rhs = [rhs_D;
       zeros(nDispDofs,1)];

[A_global, rhs, ~, ~] = assembleNeumann(face_struct, A_global, rhs, face_global_geom, lambda, mu);

%% STEP 5: SOLVE GLOBAL SYSTEM
fprintf('\n[# SOLVE THE LINEAR SYSTEM: %d x %d]\n', size(A_global,1), size(A_global,2));

switch solverOption

    case 'direct'
        sol = A_global \ rhs;
        sigma_h = sol(1:nStressDofs);
        u_h = sol(nStressDofs+1:end);

    case 'iterative'
        [sigma_h, u_h] = iterative_solver(A_global, rhs, cell_struct, face_struct);

    otherwise
        error('Unknown solver option: %s', solverOption);

end

errors = computeSolutionError(cell_struct, face_struct, sigma_h, u_h, face_global_geom, lambda, mu);

%% STEP 6: ERROR COMPUTATION / MESH VISUALIZATION

writeMeshVTU(fullfile('output', [meshName, '.vtu']), V3, cell_struct, face_struct);