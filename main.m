% close all; clear clc;

addpath('tests');
setupMRSTAuto();
G = createTwoFault();
[cell_struct, face_struct, V3, cells3D] = MRSTGridConvert(G);
cell_struct = computeVolumeQuadrature(cell_struct, face_struct, V3);

B_local = computeBE(cell_struct, face_struct, V3);

% patch/unit test for B_E
test_BE_patch(B_local, cell_struct, face_struct);
test_BE_unitcube();
