function cell_struct = computeVolumeQuadratureP(cell_struct, face_struct, V3, degree)
% cell quadrature of arbitrary polynomial degree
%
% same (cell centre, face centre, edge) tetrahedral decomposition and same
% output fields as computeVolumeQuadrature, with a degree-`degree` rule

    if nargin < 4 || isempty(degree), degree = 6; end

    [bary, wref] = quadratureRules('tet', degree);
    nq = numel(wref);

    for c = 1:numel(cell_struct)

        xE = cell_struct(c).center(:);
        cf = cell_struct(c).faces;

        quad_points = [];
        quad_weights = [];

        for lf = 1:numel(cf)

            f = cf(lf);
            xf = face_struct(f).center(:);
            verts = face_struct(f).verts(:);
            m = numel(verts);

            for k = 1:m

                vi = V3(verts(k), :).';
                vip1 = V3(verts(mod(k, m) + 1), :).';

                volumeTet = abs(det([xf - xE, vi - xE, vip1 - xE])) / 6;

                if volumeTet == 0
                    continue;
                end

                tet = [xE.'; xf.'; vi.'; vip1.'];

                quad_points = [quad_points, (bary * tet).'];      %#ok<AGROW>
                quad_weights = [quad_weights, volumeTet * wref.']; %#ok<AGROW>

            end

        end

        cell_struct(c).quad_points = quad_points;
        cell_struct(c).quad_weights = quad_weights;
        cell_struct(c).nQuad = size(quad_points, 2);
        cell_struct(c).nTet = size(quad_points, 2) / nq;

    end

end
