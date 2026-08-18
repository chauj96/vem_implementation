function face_struct = computeFaceQuadratureP(face_struct, V3, degree)
% face quadrature of arbitrary polynomial degree
%
% same fan decomposition and same output fields as computeFaceQuadrature, but
% with a degree-`degree` rule on each triangle instead of the fixed degree-2
% rule. needed because the section 5.1 Dirichlet data is degree 6, so the
% boundary integrand is degree 7

    if nargin < 3 || isempty(degree), degree = 6; end

    [bary, wref] = quadratureRules('tri', degree);
    nq = numel(wref);

    for lf = 1:numel(face_struct)

        fC = face_struct(lf).center.';
        fv = face_struct(lf).verts;
        m = numel(fv);

        quad_points = zeros(3, m*nq);
        quad_weights = zeros(1, m*nq);

        for k = 1:m

            v1 = V3(fv(k), :);
            v2 = V3(fv(mod(k, m) + 1), :);

            tri = [fC; v1; v2];

            areaTri = 0.5 * norm(cross(v1 - fC, v2 - fC));

            idx = (k-1)*nq + (1:nq);
            quad_points(:, idx) = (bary * tri).';
            quad_weights(idx) = areaTri * wref.';

        end

        face_struct(lf).quad_points = quad_points;
        face_struct(lf).quad_weights = quad_weights;
        face_struct(lf).nQuad = size(quad_points, 2);
        face_struct(lf).nTri = m;

    end

end
