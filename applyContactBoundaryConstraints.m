function [K_bc, RHS_bc, bc] = applyContactBoundaryConstraints(K_global, RHS, cell_struct, face_struct, face_global_geom)
% apply component-wise Neumann traction constraints.
%
% original system:
%   K_global * x = RHS
%
% constraints:
%   C * x = d
%
% sparse face-local elimination:
%   x = xp + Z * y
%
% reduced system:
%   K_bc * y = RHS_bc
%
% reconstruct:
%   solution = bc.xp + bc.Z * y

    nTotal = size(K_global, 1);

    assert(size(K_global, 2) == nTotal, 'K_global must be square.');
    assert(numel(RHS) == nTotal, 'RHS dimension mismatch.');

    % particular solution
    xp = zeros(nTotal, 1);

    % sparse transformation x = xp + Z * y
    % initially identity; constrained face blocks are replaced
    constrained = false(nTotal, 1);

    blockDofs = {};
    blockZ = {};

    nConstraints = 0;

    for f = 1:numel(face_struct)

        face = face_struct(f);

        if ~face.is_boundary
            continue;
        end

        mask = face.traction_mask(:);

        if ~any(mask)
            continue;
        end

        assert(numel(mask) == 3, 'Invalid traction mask on face %d.', f);

        sigmaDofs = face.stress_dofs(:);

        assert(numel(sigmaDofs) == 6, 'Expected 6 stress DOFs on face %d.', f);
        assert(~any(constrained(sigmaDofs)), 'Overlapping boundary stress DOFs on face %d.', f);

        geom = face_global_geom(f);

        xf = face.center(:);
        Af = face.area;

        qp = face.quad_points;
        w = face.quad_weights;

        % determine exterior outward normal
        e = face.cells(face.cells > 0);

        assert(numel(e) == 1, 'Boundary face %d must have one adjacent cell.', f);

        outwardDirection = xf - cell_struct(e).center(:);
        normalSign = sign(dot(geom.n, outwardDirection));

        assert(normalSign ~= 0, 'Cannot determine outward normal on face %d.', f);

        % build local traction constraint matrix C_f * sigma_f = d_f
        % enforce selected Cartesian components in the discrete L2 sense
        M = zeros(6, 6);
        g = zeros(6, 1);

        P = diag(double(mask));
        tN = face.traction_value(:);

        for q = 1:numel(w)

            Phi = evaluateFaceBasis(qp(:,q), xf, geom.Qf, geom.t1, geom.t2, geom.n, Af, geom.m20, geom.m02);

            M = M + w(q) * (Phi.' * P * Phi);
            g = g + w(q) * (Phi.' * P * tN);

        end

        % convert global-normal traction to outward traction
        Cf = normalSign * M;
        df = g;

        % remove dependent constraint equations using SVD
        [U, S, V] = svd(Cf);

        s = diag(S);

        tol = max(size(Cf)) * eps(max(s(1), 1));
        r = sum(s > tol);

        if r == 0
            continue;
        end

        % particular solution of Cf * sigma_f = df
        sigmaPart = V(:,1:r) * ((U(:,1:r).' * df) ./ s(1:r));

        % check consistency of prescribed traction
        residual = norm(Cf * sigmaPart - df);

        assert(residual < 1e-8 * max(1, norm(df)), 'Inconsistent traction constraints on face %d.', f);

        % local null-space basis
        Zf = V(:,r+1:end);

        xp(sigmaDofs) = sigmaPart;

        constrained(sigmaDofs) = true;

        blockDofs{end+1,1} = sigmaDofs;
        blockZ{end+1,1} = Zf;

        nConstraints = nConstraints + r;

    end

    % construct global sparse transformation Z
    freeDofs = find(~constrained);

    nReduced = numel(freeDofs);

    for k = 1:numel(blockZ)
        nReduced = nReduced + size(blockZ{k}, 2);
    end

    % identity mapping for unconstrained DOFs
    nFree = numel(freeDofs);

    I = freeDofs;
    J = (1:nFree).';
    V = ones(nFree, 1);

    colOffset = nFree;

    % local null-space blocks
    for k = 1:numel(blockDofs)

        dofs = blockDofs{k};
        Zf = blockZ{k};

        [ii, jj, vv] = find(Zf);

        I = [I; dofs(ii)];
        J = [J; colOffset + jj];
        V = [V; vv];

        colOffset = colOffset + size(Zf, 2);

    end

    Z = sparse(I, J, V, nTotal, nReduced);

    % reduced system
    K_bc = Z.' * K_global * Z;
    RHS_bc = Z.' * (RHS - K_global * xp);

    % reconstruction information
    bc.xp = xp;
    bc.Z = Z;
    bc.nConstraints = nConstraints;

    fprintf('\n[# NEUMANN BOUNDARY CONSTRAINTS]\n');
    fprintf('  Original DOFs    : %d\n', nTotal);
    fprintf('  Constraints      : %d\n', nConstraints);
    fprintf('  Reduced DOFs     : %d\n', nReduced);
    fprintf('  Matrix nnz       : %d\n', nnz(K_bc));

end
