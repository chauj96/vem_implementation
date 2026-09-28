function plotStressProfile(cell_struct, face_struct, sigma_h, u_h, ...
                           face_global_geom, B_local, P_local, lambda, mu)
% plot displacement and stress profiles along z
%
% self-weight column:
%   Omega = [0,1] x [0,1] x [0,2]

    nCells = numel(cell_struct);

    zCell = zeros(nCells,1);

    % numerical displacement
    ux_h = zeros(nCells,1);
    uy_h = zeros(nCells,1);
    uz_h = zeros(nCells,1);

    % numerical stress
    sxx_h = zeros(nCells,1);
    syy_h = zeros(nCells,1);
    szz_h = zeros(nCells,1);


    %% reconstruct numerical solution

    for e = 1:nCells

        xE = cell_struct(e).center(:);
        zCell(e) = xE(3);

        % -----------------------------------------
        % displacement
        % -----------------------------------------

        dDofs = 6*(e-1) + (1:6);
        uE = u_h(dDofs);

        % first three RM coefficients
        ux_h(e) = uE(1);
        uy_h(e) = uE(2);
        uz_h(e) = uE(3);


        % -----------------------------------------
        % stress
        % -----------------------------------------

        elementFaces = cell_struct(e).faces(:).';
        nElementFaces = numel(elementFaces);

        sDofs = reshape( ...
            6*(elementFaces-1) + (1:6).', [], 1);

        sigmaE_global = sigma_h(sDofs);

        % global -> element-local convention
        TE = zeros(6*nElementFaces);

        for lf = 1:nElementFaces

            f = elementFaces(lf);
            idx = 6*(lf-1) + (1:6);

            localGeom = B_local{e}.geom(lf);
            globalGeom = face_global_geom(f);

            Tf = computeFaceTransformation( ...
                f, localGeom, globalGeom, face_struct);

            TE(idx,idx) = Tf;

        end

        sigmaE_local = TE * sigmaE_global;

        % polynomial stress reconstruction
        projCoeff = P_local{e} * sigmaE_local;

        % ordering:
        % [xx, yy, zz, yz, xz, xy]
        sxx_h(e) = projCoeff(1);
        syy_h(e) = projCoeff(2);
        szz_h(e) = projCoeff(3);

    end


    %% sort cells by z

    [zCell, order] = sort(zCell);

    ux_h = ux_h(order);
    uy_h = uy_h(order);
    uz_h = uz_h(order);

    sxx_h = sxx_h(order);
    syy_h = syy_h(order);
    szz_h = szz_h(order);


    %% exact solution

    zExact = linspace(0,2,500).';

    % displacement
    ux_exact = zeros(size(zExact));
    uy_exact = zeros(size(zExact));

    uz_exact = -(4*zExact - zExact.^2) ...
               /(2*(lambda + 2*mu));

    % stress
    factor = -(2-zExact)/(lambda + 2*mu);

    sxx_exact = lambda * factor;
    syy_exact = lambda * factor;
    szz_exact = (lambda + 2*mu) * factor;


    %% plot

    figure;

    % ==========================================================
    % displacement
    % ==========================================================

    subplot(2,1,1);

    plot(zExact, ux_exact, '-', 'LineWidth', 1.5, ...
        'DisplayName', 'u_X exact');
    hold on;

    plot(zExact, uy_exact, '-', 'LineWidth', 1.5, ...
        'DisplayName', 'u_Y exact');

    plot(zExact, uz_exact, '-', 'LineWidth', 1.5, ...
        'DisplayName', 'u_Z exact');

    plot(zCell, ux_h, '-', ...
        'DisplayName', 'u_X');

    plot(zCell, uy_h, '-', ...
        'DisplayName', 'u_Y');

    plot(zCell, uz_h, '-', ...
        'DisplayName', 'u_Z');

    xlabel('$z$', 'Interpreter', 'latex');
    ylabel('$u$', 'Interpreter', 'latex');

    legend('Location', 'best');
    grid on;
    box on;


    % ==========================================================
    % stress
    % ==========================================================

    subplot(2,1,2);

    plot(zExact, sxx_exact, '-', 'LineWidth', 1.5, ...
        'DisplayName', '\sigma_{xx} exact');
    hold on;

    plot(zExact, syy_exact, '-', 'LineWidth', 1.5, ...
        'DisplayName', '\sigma_{yy} exact');

    plot(zExact, szz_exact, '-', 'LineWidth', 1.5, ...
        'DisplayName', '\sigma_{zz} exact');

    plot(zCell, sxx_h, '-', ...
        'DisplayName', '\sigma_{xx}');

    plot(zCell, syy_h, '-', ...
        'DisplayName', '\sigma_{yy}');

    plot(zCell, szz_h, '-', ...
        'DisplayName', '\sigma_{zz}');

    xlabel('$z$', 'Interpreter', 'latex');
    ylabel('$\sigma$', 'Interpreter', 'latex');

    legend('Location', 'best');
    grid on;
    box on;

end

% function plotStressProfile(cell_struct, face_struct, sigma_h, ...
%                            face_global_geom, B_local, P_local)
% % plot reconstructed sigma_zz against the exact stress profile
% %
% % self-weight column:
% %   Omega = [0,1] x [0,1] x [0,2]
% %   sigma_zz_exact(z) = z - 2
% 
%     nCells = numel(cell_struct);
% 
%     zCell = zeros(nCells,1);
%     sigmaZZ_h = zeros(nCells,1);
% 
%     for e = 1:nCells
% 
%         xE = cell_struct(e).center(:);
%         zCell(e) = xE(3);
% 
%         elementFaces = cell_struct(e).faces(:).';
%         nElementFaces = numel(elementFaces);
% 
%         % global stress DOFs associated with this element
%         sDofs = reshape( ...
%             6*(elementFaces-1) + (1:6).', [], 1);
% 
%         sigmaE_global = sigma_h(sDofs);
% 
%         % transformation:
%         % global/reference face convention -> element-local convention
%         TE = zeros(6*nElementFaces);
% 
%         for lf = 1:nElementFaces
% 
%             f = elementFaces(lf);
% 
%             idx = 6*(lf-1) + (1:6);
% 
%             localGeom = B_local{e}.geom(lf);
%             globalGeom = face_global_geom(f);
% 
%             Tf = computeFaceTransformation( ...
%                 f, localGeom, globalGeom, face_struct);
% 
%             TE(idx,idx) = Tf;
% 
%         end
% 
%         sigmaE_local = TE * sigmaE_global;
% 
%         % polynomial stress projection
%         projCoeff = P_local{e} * sigmaE_local;
% 
%         % coefficient ordering:
%         % [sigma_xx, sigma_yy, sigma_zz, sigma_yz, sigma_xz, sigma_xy]
%         sigmaZZ_h(e) = projCoeff(3);
% 
%     end
% 
%     %% exact solution
% 
%     zExact = linspace(0,2,200);
%     sigmaZZ_exact = zExact - 2;
% 
% 
%     %% plot
% 
%     figure;
% 
%     plot(zExact, sigmaZZ_exact, ...
%         'LineWidth', 1.5, ...
%         'DisplayName', 'Exact');
% 
%     hold on;
% 
%     scatter(zCell, sigmaZZ_h, ...
%         10, 'filled', ...
%         'DisplayName', 'Hybridized VEM');
% 
%     xlabel('$z$', 'Interpreter', 'latex');
%     ylabel('$\sigma_{zz}$', 'Interpreter', 'latex');
% 
%     legend('Location', 'best', ...
%            'Interpreter', 'latex');
% 
%     grid on;
%     box on;
% 
%     title('$\sigma_{zz}$ stress profile', ...
%           'Interpreter', 'latex');
% 
% 
%     sigmaZZ_exact_cell = zCell - 2;
%     sigmaZZ_error = sigmaZZ_h - sigmaZZ_exact_cell;
% 
%     figure;
% 
%     scatter(zCell, sigmaZZ_error, 12, 'filled');
%     hold on;
%     yline(0, '--', 'LineWidth', 1.2);
% 
%     xlabel('$z$', 'Interpreter', 'latex');
%     ylabel('$\sigma_{zz,h} - \sigma_{zz}^{\mathrm{exact}}$', ...
%            'Interpreter', 'latex');
% 
%     title('$\sigma_{zz}$ error', 'Interpreter', 'latex');
% 
%     grid on;
%     box on;
% 
% end