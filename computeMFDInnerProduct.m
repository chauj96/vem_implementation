function A_mfd_local = computeMFDInnerProduct(cell_struct, face_struct)
% Construct the local MFD flux mass matrix A_E
%
% A_E approximates
%
%       (K^{-1} w, z)_E
%
% using one integrated normal-flux DOF per face.

    nCells = numel(cell_struct);
    dim = 3;

    A_mfd_local = cell(nCells, 1);

    % debugging
    geom_err = zeros(nCells,1);
    sym_err  = zeros(nCells,1);
    min_eig  = zeros(nCells,1);

    for e = 1:nCells

        face_ids = cell_struct(e).faces;
        nFacesE = numel(face_ids);

        xE = cell_struct(e).center(:);
        volume = cell_struct(e).volume;
        K_perm = cell_struct(e).perm;

        C = zeros(nFacesE, dim);
        N = zeros(nFacesE, dim);
        area = zeros(nFacesE, 1);

        for lf = 1:nFacesE

            f = face_ids(lf);

            xf = face_struct(f).center(:);
            nf = face_struct(f).normal(:);
            Af = face_struct(f).area;

            s = cell_struct(e).faces_orientation(lf);

            C(lf,:) = (xf - xE).';
            N(lf,:) = (Af * s * nf).';
            area(lf) = Af;
        end

        % Simple MFD stabilization
        t = 6 * trace(K_perm) / dim;

        Q = orth(N ./ area);
        U = eye(nFacesE) - Q * Q';

        Di = diag(1 ./ area);

        A_cons = (C * (K_perm \ C')) / volume;
        A_stab = (volume / t) * (Di * U * Di);

        A_E = A_cons + A_stab;

        A_mfd_local{e} = A_E;

        % -----------------------------------------------------
        % Debugging
        % -----------------------------------------------------

        % MFD geometric consistency: C^T N = |E| I
        geom_err(e) = norm(C.'*N - volume*eye(dim), 'fro') ...
                      / max(volume, eps);

        % symmetry
        sym_err(e) = norm(A_E - A_E.', 'fro') ...
                     / max(norm(A_E, 'fro'), eps);

        % positive definiteness
        min_eig(e) = min(eig((A_E + A_E.')/2));

        % finite entries
        assert(all(isfinite(A_E(:))), ...
            'NaN/Inf detected in MFD matrix for cell %d.', e);
    end

    % ---------------------------------------------------------
    % Debug summary
    % ---------------------------------------------------------
    fprintf('\n=====================================================\n');
    fprintf('LOCAL MFD INNER PRODUCT CHECK\n');
    fprintf('=====================================================\n');
    fprintf('Number of cells             : %d\n', nCells);
    fprintf('Max geometry error          : %.3e\n', max(geom_err));
    fprintf('Max symmetry error          : %.3e\n', max(sym_err));
    fprintf('Minimum eigenvalue          : %.3e\n', min(min_eig));

    if max(geom_err) < 1e-10 && ...
       max(sym_err)  < 1e-12 && ...
       min(min_eig)  > 0

        fprintf('Status                       : PASS\n');
    else
        fprintf('Status                       : CHECK REQUIRED\n');
    end

    fprintf('=====================================================\n\n');

end

% function A_mfd_local = computeMFDInnerProduct(cell_struct, face_struct)
% % Construct the local MFD flux mass matrix A_E
% %
% % A_E approximates
% %
% %       (K^{-1} w, z)_E
% %
% % using one integrated normal-flux DOF per face.
% 
%     nCells = numel(cell_struct);
%     dim = 3;
% 
%     A_mfd_local = cell(nCells, 1);
% 
%     for e = 1:nCells
% 
%         face_ids = cell_struct(e).faces;
%         nFacesE = numel(face_ids);
% 
%         xE = cell_struct(e).center(:);
%         volume = cell_struct(e).volume;
%         K_perm = cell_struct(e).perm;
% 
%         C = zeros(nFacesE, dim);
%         N = zeros(nFacesE, dim);
%         area = zeros(nFacesE, 1);
% 
%         for lf = 1:nFacesE
% 
%             f = face_ids(lf);
% 
%             xf = face_struct(f).center(:);
%             nf = face_struct(f).normal(:);
%             Af = face_struct(f).area;
% 
%             s = cell_struct(e).faces_orientation(lf);
% 
%             % vector from cell center to face center
%             C(lf,:) = (xf - xE).';
% 
%             % outward area-weighted normal
%             N(lf,:) = (Af * s * nf).';
% 
%             area(lf) = Af;
%         end
% 
%         % Simple MFD stabilization
%         t = 6 * trace(K_perm) / dim;
% 
%         Q = orth(N ./ area);
%         U = eye(nFacesE) - Q * Q';
% 
%         Di = diag(1 ./ area);
% 
%         A_cons = (C * (K_perm \ C')) / volume;
%         A_stab = (volume / t) * (Di * U * Di);
% 
%         A_E = A_cons + A_stab;
% 
%         A_mfd_local{e} = A_E;
%     end
% end