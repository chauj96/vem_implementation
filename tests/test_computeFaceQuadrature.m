function test_computeFaceQuadrature()
% verify that computeFaceQuadrature.m function returns correct values

    fprintf('Running computeFaceQuadrature() verification...\n');
    fprintf('=====================================================\n');

    tol = 1e-14;

    % create an arbitrary triangle
    V3 = [ ...
        1.3  -0.8   2.1;
        3.7   0.4   1.5;
       -0.6   2.5   4.2];
    
    v1 = V3(1,:)';
    v2 = V3(2,:)';
    v3 = V3(3,:)';
    
    % geometry
    area = 0.5 * norm(cross(v2-v1, v3-v1));
    center = (v1 + v2 + v3)/3;
    
    face_struct_test = struct();
    face_struct_test(1).verts = [1 2 3];
    face_struct_test(1).center = center;
    face_struct_test(1).area = area;
    
    % quadrature
    face_struct_test = computeFaceQuadrature(face_struct_test, V3);
    
    qp = face_struct_test(1).quad_points;
    w = face_struct_test(1).quad_weights;
    
    %% ---------- analytic moments ----------
    % constant
    I0 = area;
    
    % first moments
    I1 = area/3 * (v1 + v2 + v3);
    
    % second moments
    s = v1 + v2 + v3;
    
    M2 = area/12 * ( ...
          s*s' ...
        + v1*v1' ...
        + v2*v2' ...
        + v3*v3');
    
    tol = 1e-14;
    
    %% constant
    Iquad = sum(w);
    err = abs(Iquad - I0);
    fprintf('constant error = %.3e\n',err);
    assert(err < tol);
    
    %% x
    Iquad = sum(w .* qp(1,:));
    err = abs(Iquad - I1(1));
    fprintf('x error = %.3e\n',err);
    assert(err < tol);
    
    %% y
    Iquad = sum(w .* qp(2,:));
    err = abs(Iquad - I1(2));
    fprintf('y error = %.3e\n',err);
    assert(err < tol);
    
    %% z
    Iquad = sum(w .* qp(3,:));
    err = abs(Iquad - I1(3));
    fprintf('z error = %.3e\n',err);
    assert(err < tol);
    
    %% x^2
    Iquad = sum(w .* qp(1,:).^2);
    err = abs(Iquad - M2(1,1));
    fprintf('x^2 error = %.3e\n',err);
    assert(err < tol);
    
    %% xy
    Iquad = sum(w .* qp(1,:) .* qp(2,:));
    err = abs(Iquad - M2(1,2));
    fprintf('xy error = %.3e\n',err);
    assert(err < tol);
    
    %% xz
    Iquad = sum(w .* qp(1,:) .* qp(3,:));
    err = abs(Iquad - M2(1,3));
    fprintf('xz error = %.3e\n',err);
    assert(err < tol);
    
    %% y^2
    Iquad = sum(w .* qp(2,:).^2);
    err = abs(Iquad - M2(2,2));
    fprintf('y^2 error = %.3e\n',err);
    assert(err < tol);
    
    %% yz
    Iquad = sum(w .* qp(2,:) .* qp(3,:));
    err = abs(Iquad - M2(2,3));
    fprintf('yz error = %.3e\n',err);
    assert(err < tol);
    
    %% z^2
    Iquad = sum(w .* qp(3,:).^2);
    err = abs(Iquad - M2(3,3));
    fprintf('z^2 error = %.3e\n',err);
    assert(err < tol);

    % % create a reference triangle
    % V3 = [ ...
    %     0 0 0;
    %     1 0 0;
    %     0 1 0];
    % 
    % face_struct_test = struct();
    % 
    % face_struct_test(1).verts = [1 2 3];
    % face_struct_test(1).center = [1/3; 1/3; 0];
    % face_struct_test(1).area = 0.5;
    % 
    % % compute quadrature
    % face_struct_test = computeFaceQuadrature(face_struct_test, V3);
    % 
    % qp = face_struct_test(1).quad_points;
    % w = face_struct_test(1).quad_weights;
    % 
    % % case 1: constant
    % Iquad = sum(w);
    % Iexact = 1/2;
    % err = abs(Iquad - Iexact);
    % fprintf('constant term error = %.3e\n', err);
    % assert(err < tol);
    % 
    % % case 2: x
    % Iquad = sum(w .* qp(1,:));
    % Iexact = 1/6;
    % err = abs(Iquad - Iexact);
    % fprintf('x term error = %.3e\n', err);
    % assert(err < tol);
    % 
    % % case 3: y
    % Iquad = sum(w .* qp(2,:));
    % Iexact = 1/6;
    % err = abs(Iquad - Iexact);
    % fprintf('y term error = %.3e\n', err);
    % assert(err < tol);
    % 
    % % case 4: x^2
    % Iquad = sum(w .* qp(1,:).^2);
    % Iexact = 1/12;
    % err = abs(Iquad - Iexact);
    % fprintf('x^2 term error = %.3e\n', err);
    % assert(err < tol);
    % 
    % % case 5: xy
    % Iquad = sum(w .* qp(1,:) .* qp(2,:));
    % Iexact = 1/24;
    % err = abs(Iquad - Iexact);
    % fprintf('xy term error = %.3e\n', err);
    % assert(err < tol);
    % 
    % % case 6: y^2
    % Iquad = sum(w .* qp(2,:).^2);
    % Iexact = 1/12;
    % err = abs(Iquad - Iexact);
    % fprintf('y^2 term error = %.3e\n', err);
    % assert(err < tol);

    fprintf('\n')
    fprintf('computeFaceQuadrature passed.\n');
    fprintf('=====================================================\n');
end