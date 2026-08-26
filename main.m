% clear; clc; close all
addpath('geometry', 'operators', 'assembly', 'solvers', 'io', 'tests', genpath('FACTORIZE'));
% setupMRSTAuto();  % for MRST automatic setup

%% STEP 1: LOAD MESH AND CHOOSE SYSTEM (mixed or hybrid), SOLVER TYPE (direct or iterative)
[meshOption, formOption, solverOption] = deal(2, 'hybrid', 'iterative');

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

%% STEP 4: (GLOBAL LEVEL) ASSEMBLE AND SOLVE GLOBAL SYSTEM

% apply Dirichlet boundary condition and body force on RHS
bodyForce = @(x) zeros(3, size(x,2));
rhs_D = assembleDirichletRHS(cell_struct, face_struct, face_global_geom);
rhs_f = -assembleForceRHS(cell_struct, bodyForce);

switch formOption

    case 'mixed'

        [sigma_h, u_h] = solveMixedVEM(cell_struct, face_struct, face_global_geom, B_local, K_cons, K_stab, rhs_D, rhs_f, lambda, mu, solverOption);

    case 'hybrid'
        [sigma_h, u_h] = solveHybridVEM(cell_struct, face_struct, face_global_geom, B_local, K_cons, K_stab, rhs_D, rhs_f, lambda, mu);

    otherwise
        error('Unknown VEM formulation: %s', formOption);

end

errors = computeSolutionError(cell_struct, face_struct, sigma_h, u_h, face_global_geom, lambda, mu);

%% STEP 6: ERROR COMPUTATION / MESH VISUALIZATION

writeMeshVTU(fullfile('output', [meshName, '.vtu']), V3, cell_struct, face_struct);