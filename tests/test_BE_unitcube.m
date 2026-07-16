function test_BE_unitcube()
% verify the implementation of the local divergence operator B_E on a unit
% cube by comparing the computed matrix with the analytical reference
% matrix derived by hand
%
% the comparison is performed entry by entry for each face block

    fprintf('\n');
    fprintf('=====================================================\n');
    fprintf('Running analytical unit cube verification...\n');
    fprintf('=====================================================\n');

    G = cartGrid([1 1 1],[1 1 1]);
    G = computeGeometry(G);

    [cell_struct,face_struct,V3,~] = MRSTGridConvert(G);

    B_local = computeBE(cell_struct,face_struct,V3);

    B_exact = analytical_BE_unitcube();

    B = B_local{1};

    tol = 1e-12;
    maxErr = 0;

    for f = 1:6

        cols = 6*(f-1)+(1:6);

        err = max(abs(B(:,cols)-B_exact(:,cols)),[],'all');

        fprintf('Face %d maximum entry error : %.3e\n',f,err);

        maxErr = max(maxErr,err);

    end

    fprintf('\n');
    fprintf('Maximum entrywise error : %.3e\n',maxErr);

    if maxErr < tol
        fprintf('Status                 : PASSED\n');
    else
        fprintf('Status                 : FAILED\n');
    end

    fprintf('=====================================================\n');

end

function B_exact = analytical_BE_unitcube()

    B_exact = zeros(6,36);
    
    % Face 1
    B_exact(:,1:6) = [
         0      0     -1      0          0          0
         1      0      0      0          0          0
         0     -1      0      0          0          0
         0      0      0   -1/6         0          0
         0    -1/2      0      0         0        1/12
       -1/2     0      0      0       1/12         0
    ];
    
    % Face 2
    B_exact(:,7:12) = [
         0      0      1      0          0          0
         1      0      0      0          0          0
         0      1      0      0          0          0
         0      0      0    1/6          0          0
         0    -1/2      0      0          0        1/12
        1/2     0      0      0      -1/12         0
    ];

    % Face 3
    B_exact(:,13:18) = [
         0     -1      0      0          0          0
         0      0     -1      0          0          0
         1      0      0      0          0          0
       -1/2     0      0      0       1/12         0
         0      0      0   -1/6         0          0
         0    -1/2      0      0         0        1/12
    ];
    
    % Face 4
    B_exact(:,19:24) = [
         0      1      0      0          0          0
         0      0      1      0          0          0
         1      0      0      0          0          0
        1/2     0      0      0      -1/12         0
         0      0      0    1/6          0          0
         0    -1/2      0      0          0        1/12
    ];
    
    % Face 5
    B_exact(:,25:30) = [
         1      0      0      0          0          0
         0     -1      0      0          0          0
         0      0     -1      0          0          0
         0    -1/2      0      0          0        1/12
       -1/2     0      0      0       1/12         0
         0      0      0   -1/6         0          0
    ];
    
    % Face 6
    B_exact(:,31:36) = [
         1      0      0      0          0          0
         0      1      0      0          0          0
         0      0      1      0          0          0
         0    -1/2      0      0          0        1/12
        1/2     0      0      0      -1/12         0
         0      0      0    1/6          0          0
    ];

end