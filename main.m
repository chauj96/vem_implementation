% clear; clc; close all
addpath('tests');
addpath(genpath('FACTORIZE'))
% setupMRSTAuto();  % for MRST automatic setup

%% LOAD MESH AND CHOOSE SOLVER TYPE (direct or iterative)
meshOption = 1;
solverOption = 'iterative'; 

if meshOption == 1

    % unit cube
    fprintf('\n[# MESH LOADED: unit cube model]\n');
    G = cartGrid([10 10 10],[1 1 1]);
    G = computeGeometry(G);
    [cell_struct, face_struct, V3, cells3D] = MRSTGridConvert(G);

elseif meshOption == 2
    % two-fault model (we need to use iterative solver!!)

    fprintf('\n[# MESH LOADED: two-fault model]\n');

    % way 1: create through MRST (takes long)
    % G = createTwoFault();
    % [cell_struct, face_struct, V3, cells3D] = MRSTGridConvert(G);

    % way 2: load preprocessed mesh geometry
    load('meshes/fault_mesh.mat'); 

    % way 3: read .vtu file
    % filename = 'meshes/fault_mesh.vtu';
    % [cell_struct, face_struct, V3, cells3D] = readVTU(filename);

    nCells = numel(cell_struct);
    nFaces = numel(face_struct);
    fprintf('3D mesh model info:\n');
    fprintf('  %d vertices\n',size(V3,1));
    fprintf('  %d cells\n',nCells);
    fprintf('  %d faces\n',nFaces);

end

%% SET UP GEOMETRY
% Lamé parameters (currently assuming homogeneous isotropic elasticity)
lambda = 1;
mu = 1;

cell_struct = assignCellProperties(cell_struct, V3, cells3D, lambda, mu);
face_struct = assignFaceProperties(face_struct); 

cell_struct = computeVolumeQuadrature(cell_struct, face_struct, V3);
face_struct = computeFaceQuadrature(face_struct, V3);

nStressDofs = 6*numel(face_struct);
nDispDofs = 6*numel(cell_struct);

%% COMPUTE VEM MATRICES (LOCAL LEVEL)
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

%% ASSEMBLE GLOBAL SYSTEM (need to do documentation and clean up)
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

%% SOLVE GLOBAL SYSTEM
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


%% MESH VISUALIZATION
if ~exist('output', 'dir')
    mkdir('output');
end

writeMeshVTU('output/mesh.vtu', V3, cell_struct, face_struct);

%% INTERNAL TESTS
% % patch/unit test for B_E
% test_BE_patch(cell_struct, face_struct, B_local);
% test_BE_unitcube();
% 
% % unit test for quadrature points and weights
% test_computeFaceQuadrature();
% test_computeVolumeQuadrature();
% 
% % unit test for projection matrix P_E
% test_projection_constant(cell_struct, face_struct, P_local, B_local);
% test_stabilization_constant(cell_struct, face_struct, B_local, K_stab);