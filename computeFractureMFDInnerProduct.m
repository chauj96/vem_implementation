function AF_mfd_local = computeFractureMFDInnerProduct(fracture)
% construct local 2D MFD flux inner product matrices.
%
% each fracture face is treated as a 2D cell.
%
% flux DOFs are integrated tangential normal fluxes across fracture edges.
%
% assumption:
%   tangential permeability K_F = I.

    nFaces = fracture.nFaces;
    dim = 2;

    AF_mfd_local = cell(nFaces, 1);

    geom_err = zeros(nFaces, 1);
    sym_err = zeros(nFaces, 1);
    min_eig = zeros(nFaces, 1);

    % Tangential coordinate system for x = constant fracture
    Q = [0 0; 1 0; 0 1];

    K_perm = eye(dim);

    for k = 1:nFaces

        edge_ids = fracture.face_edges{k};
        signs = fracture.face_edge_signs{k}(:);

        nEdgesK = numel(edge_ids);

        xF = fracture.face_centers(k,:)';
        areaF = fracture.face_areas(k);

        C = zeros(nEdgesK, dim);
        N = zeros(nEdgesK, dim);
        lengths = zeros(nEdgesK, 1);

        for le = 1:nEdgesK

            e = edge_ids(le);

            xe = fracture.edge_centers(e,:)';
            nu = fracture.edge_normals(e,:)';
            Le = fracture.edge_lengths(e);

            % Local outward co-normal
            nu_local = signs(le) * nu;

            C(le,:) = (Q' * (xe - xF))';
            N(le,:) = (Le * Q' * nu_local)';
            lengths(le) = Le;
        end

        %% Consistency matrix

        A_cons = (C * (K_perm \ C')) / areaF;

        %% Stabilization matrix

        t = 6 * trace(K_perm) / dim;

        N_unit = bsxfun(@rdivide, N, lengths);

        Qn = orth(N_unit);

        U = eye(nEdgesK) - Qn * Qn';

        Di = diag(1 ./ lengths);

        A_stab = (areaF / t) * (Di * U * Di);

        %% Local fracture MFD matrix

        A_F = A_cons + A_stab;

        AF_mfd_local{k} = A_F;

        %% Verification

        geom_err(k) = norm(C' * N - areaF * eye(dim), 'fro') / max(areaF, eps);

        sym_err(k) = norm(A_F - A_F', 'fro') / max(norm(A_F, 'fro'), eps);

        min_eig(k) = min(eig((A_F + A_F') / 2));

        assert(all(isfinite(A_F(:))), ...
            'NaN/Inf detected in fracture MFD matrix for face %d.', k);
    end

    %% Debug summary

    fprintf('\n[# LOCAL FRACTURE MFD INNER PRODUCT]\n');
    fprintf('  Number of fracture faces: %d\n', nFaces);
    fprintf('  Max geometry error: %.3e\n', max(geom_err));
    fprintf('  Max symmetry error: %.3e\n', max(sym_err));
    fprintf('  Minimum eigenvalue: %.3e\n', min(min_eig));

    if max(geom_err) < 1e-10 && max(sym_err) < 1e-12 && min(min_eig) > 0

        fprintf('  Status: PASS\n');
    else
        fprintf('  Status: CHECK REQUIRED\n');
    end
end
