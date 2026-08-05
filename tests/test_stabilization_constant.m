function test_stabilization_constant(cell_struct, face_struct, B_local, K_stab)
% verify that the stabilization vanishes for a constant stress field

    fprintf('Running stabilization consistency verification...\n');
    fprintf('=====================================================\n');

    sigma = [2 1 1;
             1 3 1;
             1 1 4];

    nCells = numel(cell_struct);
    allErr = zeros(nCells,1);

    maxErr = 0;
    worstElement = 0;
    worstResidual = [];

    tol = 1e-10;

    for e = 1:nCells

        nFaces = numel(cell_struct(e).faces);

        sigma_vem = zeros(6*nFaces,1);

        for lf = 1:nFaces

            geom = B_local{e}.geom(lf);

            f = cell_struct(e).faces(lf);
            Af = face_struct(f).area;

            n  = geom.n;
            t1 = geom.t1;
            t2 = geom.t2;

            traction = sigma*n;

            sigma_face = zeros(6,1);

            sigma_face(1) = Af * dot(traction,t1);
            sigma_face(2) = Af * dot(traction,t2);
            sigma_face(3) = Af * dot(traction,n);

            % constant stress -> rotational DOFs vanish
            sigma_face(4:6) = 0;

            cols = 6*(lf-1)+(1:6);
            sigma_vem(cols) = sigma_face;

        end

        residual = K_stab{e} * sigma_vem;

        err = norm(residual,inf);

        allErr(e) = err;

        if err > maxErr

            maxErr = err;
            worstElement = e;
            worstResidual = residual;

        end

    end

    nFail = nnz(allErr > tol);

    fprintf('Maximum residual         : %.3e\n', maxErr);
    fprintf('Cell with max error      : %d\n', worstElement);
    fprintf('Failed cells             : %d / %d\n', nFail, nCells);

    if nFail == 0
        fprintf('Status                   : PASSED\n');
    else
        fprintf('Status                   : FAILED\n');
    end

    fprintf('=====================================================');

    fprintf('\nResidual statistics\n');
    fprintf('Mean      : %.3e\n', mean(allErr));
    fprintf('Median    : %.3e\n', median(allErr));
    fprintf('90%%       : %.3e\n', prctile(allErr,90));
    fprintf('99%%       : %.3e\n', prctile(allErr,99));
    fprintf('Maximum   : %.3e\n', max(allErr));

    if nFail > 0

        fprintf('\nResidual (worst element):\n');
        disp(worstResidual)

    end

    fprintf('=====================================================\n');

end