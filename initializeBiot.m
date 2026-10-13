function [sigma_old, p_old] = initializeBiot(cell_struct, face_struct, face_global_geom)
% Initial conditions for consolidation of a soil layer
%
% Initial pore pressure:
%   p_init = alpha * varpi / (s0*(lambda+2*mu) + alpha^2)
%
% Initial total stress:
%   sigma_init = diag(-alpha*p_init, -alpha*p_init, -varpi)
%
% Coordinate convention: z increases downward.

    nCells = numel(cell_struct);
    nFaces = numel(face_struct);

    sigma_old = zeros(6*nFaces, 1);
    p_old = zeros(nCells, 1);

    % initialize pressure

    varpi = 1.0;

    lambda = cell_struct(1).lambda;
    mu = cell_struct(1).mu;
    alpha = cell_struct(1).alpha;
    s0 = cell_struct(1).s0;

    H = lambda + 2*mu;

    p_init = alpha * varpi / (s0*H + alpha^2);

    p_old(:) = p_init;

    % initialize total stress

    eps_zz = (-varpi + alpha*p_init)/H;
    sigma_lateral = lambda*eps_zz - alpha*p_init;
    
    sigma_init = diag([sigma_lateral, sigma_lateral, -varpi ]);

    for f = 1:nFaces

        geom = face_global_geom(f);

        % initial traction in the global face orientation
        traction_init = sigma_init * geom.n;

        % constant traction coefficients in the six-mode basis
        xf = face_struct(f).center(:);
        Af = face_struct(f).area;

        Phi = evaluateFaceBasis(xf, xf, geom.Qf, geom.t1, geom.t2, geom.n, Af, geom.m20, geom.m02);

        % constant traction is represented by the first three modes of the face traction basis.
        sigma_face = zeros(6,1);
        sigma_face(1:3) = Phi(:,1:3) \ traction_init;

        globalIdx = 6*(f-1) + (1:6);

        sigma_old(globalIdx) = sigma_face;

    end

end
