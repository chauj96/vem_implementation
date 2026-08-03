function test_BE_patch(cell_struct, face_struct, B_local)
% verify the local divergence operator B_E using consistency patch tests
%
% test 1: sigma = I
%
% test 2: sigma = [0 1 0; 1 0 0; 0 0 0]
%
% both exact divergence fields are identically zero

    tol = 1e-12;

    fprintf('\n');
    fprintf('=====================================================\n');
    fprintf('Running B_E consistency patch tests...\n');
    fprintf('=====================================================\n');

    runPatchTest(B_local, cell_struct, face_struct, ...
        eye(3), 'Identity tensor', tol);

    sigma = [0 1 0;
             1 0 0;
             0 0 0];

    runPatchTest(B_local, cell_struct, face_struct, ...
        sigma, 'Constant symmetric tensor', tol);

end


function runPatchTest(B_local, cell_struct, face_struct, ...
                      sigma, testName, tol)

    nCells = numel(B_local);

    failedCells = [];
    maxResidual = 0;
    worstCell = 0;

    for e = 1:nCells

        % matrix and element-local face geometry used by computeBE
        B_E = B_local{e}.matrix;
        geom = B_local{e}.geom;

        faces = cell_struct(e).faces;
        nFaces = numel(faces);

        dof = zeros(6*nFaces, 1);

        for lf = 1:nFaces

            f = faces(lf);
            Af = face_struct(f).area;

            % use exactly the same local frame as computeBE
            t1 = geom(lf).t1;
            t2 = geom(lf).t2;
            n  = geom(lf).n;

            traction = sigma * n;

            alpha = dot(traction, t1);
            beta  = dot(traction, t2);
            gamma = dot(traction, n);

            cols = 6*(lf-1) + (1:6);

            dof(cols) = [
                Af * alpha
                Af * beta
                Af * gamma
                0
                0
                0
            ];

        end

        divCoeff = B_E * dof;
        residual = norm(divCoeff);

        if residual > maxResidual
            maxResidual = residual;
            worstCell = e;
        end

        if residual > tol
            failedCells(end+1) = e;
        end

    end

    fprintf('\n%s\n', testName);
    fprintf('-----------------------------------------------\n');
    fprintf('Maximum residual      : %.3e\n', maxResidual);
    fprintf('Cell with max error   : %d\n', worstCell);
    fprintf('Failed cells          : %d / %d\n', ...
        numel(failedCells), nCells);

    if isempty(failedCells)
        fprintf('Status                : PASSED\n');
    else
        fprintf('Status                : FAILED\n');
    end

end