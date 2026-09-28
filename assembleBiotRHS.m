function [f_sigma, f_u, f_w, f_p] = assembleBiotRHS(cell_struct, face_struct, face_global_geom, t)
% assemble RHS vectors for Biot temporal convergence test
%
% Exact solution:
%   u = exp(t) * [x^2; y^2; z^2]
%   p = exp(t) * (x^2 + y^2 + z^2)
%
% Exact stress:
%   sigma = exp(t) * sigma_0
%
% Exact Darcy flux:
%   w = exp(t) * [-2*x; -2*y; -2*z]
%
% Source terms:
%   b = exp(t) * [2*x-6; 2*y-6; 2*z-6]
%   psi = exp(t) * (x^2+y^2+z^2 + 2*x+2*y+2*z - 6)
%
% Boundary conditions:
%   u = u_exact(x,t) on entire boundary
%   p = p_exact(x,t) on entire boundary

    et = exp(t);

    nCells = numel(cell_struct);
    nFaces = numel(face_struct);

    nStressDofs = 6 * nFaces;
    nDispDofs = 6 * nCells;
    nFluxDofs = nFaces;
    nPressureDofs = nCells;

    f_sigma = zeros(nStressDofs, 1);
    f_u = zeros(nDispDofs, 1);
    f_w = zeros(nFluxDofs, 1);
    f_p = zeros(nPressureDofs, 1);


    %% 1. BODY FORCE
    % b = exp(t) * [2*x-6; 2*y-6; 2*z-6]
    %
    % Since b is affine and xE is the cell centroid,
    % its translational moments are exactly volume*b(xE).
    % The rotational moments vanish for this particular b.

    for e = 1:nCells

        volume = cell_struct(e).volume;
        xE = cell_struct(e).center(:);

        bE = et * [ ...
            2*xE(1) - 6;
            2*xE(2) - 6;
            2*xE(3) - 6 ];

        f_u_local = volume * [bE; 0; 0; 0];

        dispDofs = 6*(e-1) + (1:6);

        f_u(dispDofs) = f_u(dispDofs) + f_u_local;

    end


    %% 2. DISPLACEMENT BOUNDARY TERM
    % f_sigma_j = integral_f g_u . (phi_j n_out) dS

    % find adjacent cell for each displacement boundary face
    faceCell_u = zeros(nFaces, 1);

    for e = 1:nCells

        elementFaces = cell_struct(e).faces;

        for lf = 1:numel(elementFaces)

            f = elementFaces(lf);

            if face_struct(f).is_displacement_bc
                faceCell_u(f) = e;
            end

        end
    end


    for f = 1:nFaces

        if ~face_struct(f).is_displacement_bc
            continue;
        end

        xf = face_struct(f).center(:);
        Af = face_struct(f).area;

        geom = face_global_geom(f);

        qp = face_struct(f).quad_points;
        qw = face_struct(f).quad_weights;

        % physical outward orientation
        e = faceCell_u(f);

        assert(e ~= 0, ...
            'Could not find adjacent cell for displacement face %d.', f);

        xE = cell_struct(e).center(:);

        outwardDirection = xf - xE;

        normalSign = sign(dot(geom.n, outwardDirection));

        assert(normalSign ~= 0, ...
            'Could not determine outward orientation for face %d.', f);

        f_sigma_face = zeros(6, 1);

        for q = 1:numel(qw)

            xq = qp(:,q);

            % exact displacement
            gu = et * [xq(1)^2; xq(2)^2; xq(3)^2];

            % traction basis in global/reference convention
            PhiG = evaluateFaceBasis( ...
                xq, xf, geom.Qf, geom.t1, geom.t2, geom.n, ...
                Af, geom.m20, geom.m02);

            % convert to physical outward convention
            PhiOut = normalSign * PhiG;

            f_sigma_face = f_sigma_face ...
                         + qw(q) * (PhiOut' * gu);

        end

        stressDofs = 6*(f-1) + (1:6);

        f_sigma(stressDofs) = ...
            f_sigma(stressDofs) + f_sigma_face;

    end


    %% 3. PRESSURE BOUNDARY TERM
    % p = x^2 + y^2 + z^2

    % find the adjacent cell for each boundary face
    faceCell = zeros(nFaces, 1);

    for e = 1:nCells

        elementFaces = cell_struct(e).faces;

        for lf = 1:numel(elementFaces)

            f = elementFaces(lf);

            if face_struct(f).is_pressure_bc
                faceCell(f) = e;
            end

        end
    end


    for f = 1:nFaces

        if ~face_struct(f).is_pressure_bc
            continue;
        end

        qp = face_struct(f).quad_points;
        qw = face_struct(f).quad_weights;
        Af = face_struct(f).area;

        % global/reference normal
        geom = face_global_geom(f);

        % physical outward orientation
        e = faceCell(f);

        assert(e ~= 0, ...
            'Could not find adjacent cell for pressure face %d.', f);

        xf = face_struct(f).center(:);
        xE = cell_struct(e).center(:);

        outwardDirection = xf - xE;

        normalSign = sign(dot(geom.n, outwardDirection));

        assert(normalSign ~= 0, ...
            'Could not determine outward orientation for face %d.', f);

        % integrate prescribed pressure
        gp_integral = 0;

        for q = 1:numel(qw)

            xq = qp(:,q);

            x = xq(1);
            y = xq(2);
            z = xq(3);

            gp = et * (x^2 + y^2 + z^2);

            gp_integral = gp_integral + qw(q) * gp;

        end

        % convert physical-outward boundary term
        % to global/reference flux convention
        f_w(f) = -normalSign * gp_integral / Af;

    end


    %% 4. FLUID SOURCE
    for e = 1:nCells
    
        qp = cell_struct(e).quad_points;
        qw = cell_struct(e).quad_weights;
    
        fp_local = 0;
    
        for q = 1:numel(qw)
    
            xq = qp(:,q);
    
            x = xq(1);
            y = xq(2);
            z = xq(3);
    
            psi = et * ( ...
                x^2 + y^2 + z^2 ...
                + 2*x + 2*y + 2*z ...
                - 6 );
    
            fp_local = fp_local + qw(q) * psi;
    
        end
    
        f_p(e) = fp_local;
    
    end

end

% function [f_sigma, f_u, f_w, f_p] = assembleBiotRHS(cell_struct, face_struct, face_global_geom)
% % assemble RHS vectors for Biot spatial convergence test
% %
% % Exact solution:
% %   u = [x^2; y^2; z^2]
% %   p = x + y + z
% %
% % Exact stress:
% %   sigma = [5*x+y+z,       0,       0;
% %                   0, x+5*y+z,       0;
% %                   0,       0, x+y+5*z]
% %
% % Exact Darcy flux:
% %   w = [-1; -1; -1]
% %
% % Source terms:
% %   b   = [-5; -5; -5]
% %   psi = 0
% %
% % Boundary conditions:
% %   u = u_exact on entire boundary
% %   p = p_exact on entire boundary
% 
%     nCells = numel(cell_struct);
%     nFaces = numel(face_struct);
% 
%     nStressDofs = 6 * nFaces;
%     nDispDofs = 6 * nCells;
%     nFluxDofs = nFaces;
%     nPressureDofs = nCells;
% 
%     f_sigma = zeros(nStressDofs, 1);
%     f_u = zeros(nDispDofs, 1);
%     f_w = zeros(nFluxDofs, 1);
%     f_p = zeros(nPressureDofs, 1);
% 
% 
%     %% 1. BODY FORCE
%     % f_u = integral_E b . R_i dE
% 
%     b = [-5; -5; -5];
% 
%     for e = 1:nCells
% 
%         volume = cell_struct(e).volume;
% 
%         f_u_local = volume * [b; 0; 0; 0];
% 
%         dispDofs = 6*(e-1) + (1:6);
% 
%         f_u(dispDofs) = f_u(dispDofs) + f_u_local;
% 
%     end
% 
% 
%     %% 2. DISPLACEMENT BOUNDARY TERM
%     % f_sigma_j = integral_f g_u . (phi_j n_out) dS
% 
%     % find adjacent cell for each displacement boundary face
%     faceCell_u = zeros(nFaces, 1);
% 
%     for e = 1:nCells
% 
%         elementFaces = cell_struct(e).faces;
% 
%         for lf = 1:numel(elementFaces)
% 
%             f = elementFaces(lf);
% 
%             if face_struct(f).is_displacement_bc
%                 faceCell_u(f) = e;
%             end
% 
%         end
%     end
% 
% 
%     for f = 1:nFaces
% 
%         if ~face_struct(f).is_displacement_bc
%             continue;
%         end
% 
%         xf = face_struct(f).center(:);
%         Af = face_struct(f).area;
% 
%         geom = face_global_geom(f);
% 
%         qp = face_struct(f).quad_points;
%         qw = face_struct(f).quad_weights;
% 
%         % physical outward orientation
%         e = faceCell_u(f);
% 
%         assert(e ~= 0, 'Could not find adjacent cell for displacement face %d.', f);
% 
%         xE = cell_struct(e).center(:);
% 
%         outwardDirection = xf - xE;
% 
%         normalSign = sign(dot(geom.n, outwardDirection));
% 
%         assert(normalSign ~= 0, 'Could not determine outward orientation for face %d.', f);
% 
%         f_sigma_face = zeros(6, 1);
% 
%         for q = 1:numel(qw)
% 
%             xq = qp(:,q);
% 
%             % exact displacement
%             gu = [xq(1)^2; xq(2)^2; xq(3)^2];
% 
%             % traction basis in global/reference convention
%             PhiG = evaluateFaceBasis(xq, xf, geom.Qf, geom.t1, geom.t2, geom.n, Af, geom.m20, geom.m02);
% 
%             % convert to physical outward convention
%             PhiOut = normalSign * PhiG;
% 
%             f_sigma_face = f_sigma_face + qw(q) * (PhiOut' * gu);
% 
%         end
% 
%         stressDofs = 6*(f-1) + (1:6);
% 
%         f_sigma(stressDofs) = f_sigma(stressDofs) + f_sigma_face;
% 
%     end
% 
% 
%     %% 3. PRESSURE BOUNDARY TERM
%     % Pressure boundary contribution in the global/reference
%     % face-normal convention.
% 
%     % find the adjacent cell for each boundary face
%     faceCell = zeros(nFaces, 1);
% 
%     for e = 1:nCells
% 
%         elementFaces = cell_struct(e).faces;
% 
%         for lf = 1:numel(elementFaces)
% 
%             f = elementFaces(lf);
% 
%             if face_struct(f).is_pressure_bc
%                 faceCell(f) = e;
%             end
% 
%         end
%     end
% 
% 
%     for f = 1:nFaces
% 
%         if ~face_struct(f).is_pressure_bc
%             continue;
%         end
% 
%         qp = face_struct(f).quad_points;
%         qw = face_struct(f).quad_weights;
%         Af = face_struct(f).area;
% 
%         % global/reference normal
%         geom = face_global_geom(f);
% 
%         % physical outward orientation
%         e = faceCell(f);
% 
%         assert(e ~= 0, 'Could not find adjacent cell for pressure face %d.', f);
% 
%         xf = face_struct(f).center(:);
%         xE = cell_struct(e).center(:);
% 
%         outwardDirection = xf - xE;
% 
%         normalSign = sign(dot(geom.n, outwardDirection));
% 
%         assert(normalSign ~= 0, 'Could not determine outward orientation for face %d.', f);
% 
%         % integrate prescribed pressure
%         gp_integral = 0;
% 
%         for q = 1:numel(qw)
% 
%             xq = qp(:,q);
% 
%             x = xq(1);
%             y = xq(2);
%             z = xq(3);
%             gp = x + y + z;
% 
%             gp_integral = gp_integral + qw(q) * gp;
% 
%         end
% 
%         % convert physical-outward boundary term to global/reference flux convention
%         f_w(f) = -normalSign * gp_integral / Af;
% 
%     end
% 
%     %% 4. FLUID SOURCE
%     % psi = 0
%     %
%     % Therefore f_p remains zero.
% 
% end