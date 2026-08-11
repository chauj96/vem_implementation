% clear; clc; close all
addpath('tests');
% setupMRSTAuto();  % for MRST automatic setup

%% LOAD MESH
meshOption = 2;

if meshOption == 1

    % unit cube
    G = cartGrid([3 3 3],[1 1 1]);
    G = computeGeometry(G);
    [cell_struct, face_struct, V3, cells3D] = MRSTGridConvert(G);

elseif meshOption == 2
    % two-fault model

    % create through MRST (takes long)
    % G = createTwoFault();
    % [cell_struct, face_struct, V3, cells3D] = MRSTGridConvert(G);

    % load preprocessed mesh geometry
    load('meshes/fault_mesh.mat');    
    fprintf('\nMesh loaded: two-fault model\n');

    % read .vtu file
    % filename = 'meshes/fault_mesh.vtu';
    % readVTU(filename);

end

%% SET UP GEOMETRY
% Lamé parameters (currently assuming homogeneous isotropic elasticity)
lambda = 1;
mu = 1;

% [cell_struct, face_struct, V3, cells3D] = MRSTGridConvert(G);
cell_struct = assignCellProperties(cell_struct, V3, cells3D, lambda, mu);
cell_struct = computeVolumeQuadrature(cell_struct, face_struct, V3);
face_struct = computeFaceQuadrature(face_struct, V3);

%% COMPUTE VEM MATRICES (LOCAL LEVEL)

% INERTIA MATRIX (quadrature approach, Tonon approach)
% A_local = computeAEQuadrature(cell_struct);
A_local = computeAE_Tonon(cell_struct, face_struct, V3);

% DIVERGENCE MATRIX (quadrature approach, closed form approach from notes)
% B_local = computeBEQuadrature(cell_struct, face_struct, V3);
B_local = computeBE(cell_struct, face_struct, V3);

% PROJECTION MATRIX
D_local = computeDE(cell_struct, A_local, B_local);
Bproj = computeBproj(cell_struct, face_struct, D_local, B_local);
P_local = computeProjectionMatrix(cell_struct, Bproj);

K_cons = computeConsistency(cell_struct, P_local);
K_stab = computeStabilization(cell_struct, face_struct, B_local, P_local);

%% MESH VISUALIZATION
if ~exist('output', 'dir')
    mkdir('output');
end

writeMeshVTU('output/twoFault_mesh.vtu', V3, cell_struct, face_struct);

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