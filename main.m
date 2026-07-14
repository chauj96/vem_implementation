close all; clear clc;

setupMRSTAuto();
G = createTwoFault();
[cell_struct, face_struct, V3, cells3D] = MRSTGridConvert(G);

B_local = computeBE(cell_struct, face_struct, V3);

e = 1;

B_E = B_local{e};

fprintf('Element %d\n', e);
fprintf('Number of faces: %d\n', ...
    numel(cell_struct(e).faces));

fprintf('Size of B_E: %d x %d\n', ...
    size(B_E,1), size(B_E,2));

disp(B_E);