close all; clear clc;
addpath('tests');
% setupMRSTAuto();  % for MRST automatic setup

% (option 1) unit cube for simplicity
% G = cartGrid([10 10 10],[1 1 1]);
% G = computeGeometry(G);

% (option 2) two fault model
G = createTwoFault(); 

[cell_struct, face_struct, V3, cells3D] = MRSTGridConvert(G);
cell_struct = computeVolumeQuadrature(cell_struct, face_struct, V3);
face_struct = computeFaceQuadrature(face_struct, V3);

% A_quad = computeAEQuadrature(cell_struct);  
A_local = computeAE_Tonon(cell_struct, face_struct, V3);

% B_quad = computeBEQuadrature(cell_struct, face_struct, V3);
B_local = computeBE(cell_struct, face_struct, V3);

D_local = computeDE(cell_struct, A_local, B_quad);
Bproj = computeBproj(cell_struct, face_struct, D_local, B_local);
P_local = computeProjectionMatrix(cell_struct, Bproj);

%% TESTING
% patch/unit test for B_E
test_BE_patch(cell_struct, face_struct, B_local);
test_BE_unitcube();

% unit test for quadrature points and weights
test_computeFaceQuadrature();
test_computeVolumeQuadrature();

% unit test for projection matrix P_E
test_projection_constant(cell_struct, face_struct, P_local, B_local);