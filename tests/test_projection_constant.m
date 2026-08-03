function test_projection_constant(cell_struct, face_struct, P_local, B_local)
% verify the polynomial consistency of the local projection operator

    fprintf('\n');
    fprintf('Running projection consistency verification...\n');
    fprintf('=====================================================\n');

    sigma = [2 1 1;
             1 3 1;
             1 1 4];

    sigma_sym = [2;
                 3;
                 4;
                 sqrt(2);
                 sqrt(2);
                 sqrt(2)];

    nCells = numel(cell_struct);

    allErr = zeros(nCells,1);

    maxErr = 0;
    worstElement = 0;
    worstCoeff = [];

    tol = 1e-10;

    for e = 1:nCells

        nFaces = numel(cell_struct(e).faces);

        sigma_vem = zeros(6 * nFaces,1);

        for lf = 1:nFaces

            geom = B_local{e}.geom(lf);

            f = cell_struct(e).faces(lf);
            Af = face_struct(f).area;

            n = geom.n;
            t1 = geom.t1;
            t2 = geom.t2;

            traction = sigma*n;

            sigma_face = zeros(6,1);

            sigma_face(1) = Af * dot(traction,t1);
            sigma_face(2) = Af * dot(traction,t2);
            sigma_face(3) = Af * dot(traction,n);

            % constant stress -> rotational DOFs vanish
            sigma_face(4:6) = 0;

            cols = 6 * (lf-1) + (1:6);
            sigma_vem(cols) = sigma_face;

        end

        coeff = P_local{e} * sigma_vem;

        err = norm(coeff - sigma_sym, inf);

        allErr(e) = err;

        if err > maxErr
            maxErr = err;
            worstElement = e;
            worstCoeff = coeff;
        end

    end

    nFail = nnz(allErr > tol);

    fprintf('Maximum coefficient error : %.3e\n', maxErr);
    fprintf('Cell with max error       : %d\n', worstElement);
    fprintf('Failed cells              : %d / %d\n', nFail, nCells);

    if nFail == 0
        fprintf('Status                    : PASSED\n');
    else
        fprintf('Status                    : FAILED\n');
    end

    fprintf('=====================================================\n');

    fprintf('\nError statistics\n');
    fprintf('Mean      : %.3e\n', mean(allErr));
    fprintf('Median    : %.3e\n', median(allErr));
    fprintf('90%%       : %.3e\n', prctile(allErr,90));
    fprintf('99%%       : %.3e\n', prctile(allErr,99));
    fprintf('Maximum   : %.3e\n', max(allErr));

    if nFail > 0

        fprintf('\nRecovered coefficients (worst element):\n');
        disp(worstCoeff)

        fprintf('Expected coefficients:\n');
        disp(sigma_sym)

        fprintf('Difference:\n');
        disp(worstCoeff - sigma_sym)

    end

    fprintf('=====================================================\n');

end