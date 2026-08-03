function test_computeVolumeQuadrature()
% verify that computeVolumeQuadrature.m integrates quadratic polynomials
% exactly on the unit cube.

    fprintf('Running computeVolumeQuadrature() verification...\n');
    fprintf('=====================================================\n');

    tol = 1e-14;

    % create unit cube

    % setupMRSTAuto();

    G = cartGrid([1 1 1],[1 1 1]);
    G = computeGeometry(G);

    [cell_struct, face_struct, V3, ~] = MRSTGridConvert(G);

    % compute quadrature
    cell_struct = computeVolumeQuadrature(cell_struct, face_struct, V3);

    qp = cell_struct(1).quad_points;
    w = cell_struct(1).quad_weights;

    % exact integrals over unit cube
    I0 = 1.0;

    Ix = 1/2;
    Iy = 1/2;
    Iz = 1/2;

    Ixx = 1/3;
    Ixy = 1/4;
    Ixz = 1/4;

    Iyy = 1/3;
    Iyz = 1/4;

    Izz = 1/3;

    % constant
    Iquad = sum(w);
    err = abs(Iquad - I0);
    fprintf('constant error = %.3e\n', err);
    assert(err < tol);

    % x
    Iquad = sum(w .* qp(1,:));
    err = abs(Iquad - Ix);
    fprintf('x error = %.3e\n', err);
    assert(err < tol);

    % y
    Iquad = sum(w .* qp(2,:));
    err = abs(Iquad - Iy);
    fprintf('y error = %.3e\n', err);
    assert(err < tol);

    % z
    Iquad = sum(w .* qp(3,:));
    err = abs(Iquad - Iz);
    fprintf('z error = %.3e\n', err);
    assert(err < tol);

    % x^2
    Iquad = sum(w .* qp(1,:).^2);
    err = abs(Iquad - Ixx);
    fprintf('x^2 error = %.3e\n', err);
    assert(err < tol);

    % xy
    Iquad = sum(w .* qp(1,:) .* qp(2,:));
    err = abs(Iquad - Ixy);
    fprintf('xy error = %.3e\n', err);
    assert(err < tol);

    % xz
    Iquad = sum(w .* qp(1,:) .* qp(3,:));
    err = abs(Iquad - Ixz);
    fprintf('xz error = %.3e\n', err);
    assert(err < tol);

    % y^2
    Iquad = sum(w .* qp(2,:).^2);
    err = abs(Iquad - Iyy);
    fprintf('y^2 error = %.3e\n', err);
    assert(err < tol);

    % yz
    Iquad = sum(w .* qp(2,:) .* qp(3,:));
    err = abs(Iquad - Iyz);
    fprintf('yz error = %.3e\n', err);
    assert(err < tol);

    % z^2
    Iquad = sum(w .* qp(3,:).^2);
    err = abs(Iquad - Izz);
    fprintf('z^2 error = %.3e\n', err);
    assert(err < tol);

    fprintf('\n');
    fprintf('computeVolumeQuadrature passed.\n');
    fprintf('=====================================================\n');

end