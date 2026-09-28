function [sigma_old, p_old] = initializeBiot(cell_struct, face_struct, face_global_geom)
% initialize exact discrete stress and pressure for Biot convergence test
%
% Exact solution:
%   u = [x^2; y^2; z^2]
%   p = x^2 + y^2 + z^2

    nCells = numel(cell_struct);
    nFaces = numel(face_struct);

    sigma_old = zeros(6*nFaces, 1);
    p_old = zeros(nCells, 1);


    %% 1. EXACT PRESSURE DOFs
    %
    % Pressure DOF = cell average:
    %
    %   p_E = (1/|E|) * integral_E p dE

    for e = 1:nCells

        qp = cell_struct(e).quad_points;
        qw = cell_struct(e).quad_weights;
        volume = cell_struct(e).volume;

        p_integral = 0;

        for q = 1:numel(qw)

            xq = qp(:,q);

            x = xq(1);
            y = xq(2);
            z = xq(3);

            p_exact = x^2 + y^2 + z^2;

            p_integral = p_integral + qw(q) * p_exact;

        end

        p_old(e) = p_integral / volume;

    end


    %% 2. EXACT STRESS DOFs
    %
    % On each face:
    %
    %   t_exact(x) = sigma_exact(x) * n_global
    %
    % Represent this traction in the six global face traction modes.

    for f = 1:nFaces

        xf = face_struct(f).center(:);
        Af = face_struct(f).area;

        geom = face_global_geom(f);

        qp = face_struct(f).quad_points;
        qw = face_struct(f).quad_weights;

        % L2 Gram matrix and RHS
        Gface = zeros(6,6);
        bface = zeros(6,1);

        for q = 1:numel(qw)

            xq = qp(:,q);

            x = xq(1);
            y = xq(2);
            z = xq(3);

            % exact pressure
            p = x^2 + y^2 + z^2;

            % exact stress
            sigma_exact = [ ...
                6*x + 2*y + 2*z - p,  0,                      0;
                0,                      2*x + 6*y + 2*z - p,  0;
                0,                      0,                      2*x + 2*y + 6*z - p ];

            % exact traction using global/reference normal
            traction_exact = sigma_exact * geom.n;

            % global/reference traction basis
            Phi = evaluateFaceBasis( ...
                xq, xf, geom.Qf, geom.t1, geom.t2, geom.n, ...
                Af, geom.m20, geom.m02);

            % L2 projection onto face traction space
            Gface = Gface + qw(q) * (Phi' * Phi);

            bface = bface + qw(q) * (Phi' * traction_exact);

        end

        sigma_face = Gface \ bface;

        globalIdx = 6*(f-1) + (1:6);

        sigma_old(globalIdx) = sigma_face;

    end

end

% function [sigma_old, p_old] = initializeBiot(cell_struct, face_struct, face_global_geom)
% % initialize exact discrete stress and pressure for Biot convergence test
% 
%     nCells = numel(cell_struct);
%     nFaces = numel(face_struct);
% 
%     sigma_old = zeros(6*nFaces, 1);
%     p_old = zeros(nCells, 1);
% 
% 
%     %% 1. EXACT PRESSURE DOFs
%     for e = 1:nCells
% 
%         xE = cell_struct(e).center(:);
% 
%         p_old(e) = xE(1) + xE(2) + xE(3);
% 
%     end
% 
% 
%     %% 2. EXACT STRESS DOFs
%     %
%     % On each face:
%     %
%     %   t_exact(x) = sigma_exact(x) * n_global
%     %
%     % Represent this traction in the six global face traction modes.
% 
%     for f = 1:nFaces
% 
%         xf = face_struct(f).center(:);
%         Af = face_struct(f).area;
% 
%         geom = face_global_geom(f);
% 
%         qp = face_struct(f).quad_points;
%         qw = face_struct(f).quad_weights;
% 
%         % L2 Gram matrix and RHS
%         Gface = zeros(6,6);
%         bface = zeros(6,1);
% 
%         for q = 1:numel(qw)
% 
%             xq = qp(:,q);
% 
%             % exact stress
%             x = xq(1);
%             y = xq(2);
%             z = xq(3);
% 
%             sigma_exact = [ ...
%                             5*x + y + z,  0,              0;
%                             0,              x + 5*y + z,  0;
%                             0,              0,              x + y + 5*z ];
% 
%             % exact traction using global/reference normal
%             traction_exact = sigma_exact * geom.n;
% 
%             % global/reference traction basis
%             Phi = evaluateFaceBasis(xq, xf, geom.Qf, geom.t1, geom.t2, geom.n, Af, geom.m20, geom.m02);
% 
%             % L2 projection onto face traction space
%             Gface = Gface + qw(q) * (Phi' * Phi);
% 
%             bface = bface + qw(q) * (Phi' * traction_exact);
% 
%         end
% 
%         sigma_face = Gface \ bface;
% 
%         globalIdx = 6*(f-1) + (1:6);
% 
%         sigma_old(globalIdx) = sigma_face;
% 
%     end
% 
% end