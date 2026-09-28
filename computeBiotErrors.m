function err = computeBiotErrors(cell_struct, face_struct, face_global_geom, V3, B_local, D_local, P_local, A_mfd_local, sigma_h, u_h, w_h, p_h, problem)
% Spatial error norms for the fully mixed Biot discretization
%
% Mechanics:
%   E_u      = || u - u_h ||_0
%   E_sigdiv = || div(sigma) - div(sigma_h) ||_0
%   E_sigPi  = || sigma - Pi_E sigma_h ||_0
%   E_sig    = traction-based discrete stress error
%
% Flow:
%   E_p      = || p - p_h ||_0
%   E_w      = discrete MFD flux error induced by A_E
%
% All errors are reported as relative errors whenever the exact
% reference norm is nonzero.


    nCells = numel(cell_struct);
    nFaces = numel(face_struct);

    eU2   = 0;   nU2   = 0;
    eDiv2 = 0;   nDiv2 = 0;
    ePi2  = 0;   nPi2  = 0;
    eP2   = 0;   nP2   = 0;
    eW2   = 0;   nW2   = 0;

    s2 = 1/sqrt(2);

    % CELL ERRORS
    for e = 1:nCells

        qp = cell_struct(e).quad_points;
        qw = cell_struct(e).quad_weights(:).';

        xE = cell_struct(e).center(:);
        r = qp - xE;

        % Local stress DOFs in element convention
        sl = localStressDofs(e, cell_struct, face_struct, face_global_geom, B_local, sigma_h);

        % 1. DISPLACEMENT ERROR
        % u_h |_E = a + omega x (x - x_E)

        a = u_h(6*(e-1) + (1:3));
        omega = u_h(6*(e-1) + (4:6));

        uh = a + cross(repmat(omega, 1, size(r,2)), r, 1);
        ue = problem.u(qp);

        eU2 = eU2 + sum(sum((ue - uh).^2, 1) .* qw);
        nU2 = nU2 + sum(sum(ue.^2, 1) .* qw);

        % 2. STRESS DIVERGENCE ERROR
        % div(sigma_h)|_E = alpha_E + omega_E x (x - x_E)
        d = D_local{e} * sl;

        dh = d(1:3) + cross(repmat(d(4:6), 1, size(r,2)), r, 1);

        de = problem.divsigma(qp);

        eDiv2 = eDiv2 + sum(sum((de - dh).^2, 1) .* qw);
        nDiv2 = nDiv2 + sum(sum(de.^2, 1) .* qw);

        % 3. PROJECTED STRESS ERROR
        % Pi_E sigma_h is constant on each cell
        c = P_local{e} * sl;

        Sh = [     c(1), s2*c(6), s2*c(5); ...
                s2*c(6),    c(2), s2*c(4); ...
                s2*c(5), s2*c(4),    c(3)];

        Se = problem.sigma(qp);

        dS = Se - Sh;

        ePi2 = ePi2 + sum(reshape(sum(sum(dS.^2, 1), 2), 1, []) .* qw);

        nPi2 = nPi2 + sum(reshape(sum(sum(Se.^2, 1), 2), 1, []) .* qw);

        % 4. PRESSURE ERROR
        % p_h is one cellwise constant DOF per element
        ph = p_h(e);
        pe = problem.p(qp);

        pe = reshape(pe, 1, []);

        eP2 = eP2 + sum((pe - ph).^2 .* qw);
        nP2 = nP2 + sum(pe.^2 .* qw);

        % 5. MFD FLUX ERROR
        % Local MFD DOF:
        %
        %   w_E,f = integral_f w . n_E,f df
        face_ids = cell_struct(e).faces(:);
        orient = cell_struct(e).faces_orientation(:);

        % numerical local flux DOFs
        wh_local = orient .* w_h(face_ids);

        % exact interpolated local flux DOFs
        wI_local = zeros(numel(face_ids), 1);

        for lf = 1:numel(face_ids)

            f = face_ids(lf);

            qp_f = face_struct(f).quad_points;
            qw_f = face_struct(f).quad_weights(:).';

            nf = face_struct(f).normal(:);

            % cell-outward normal
            nout = orient(lf) * nf;

            % exact Darcy flux at face quadrature points
            we = problem.w(qp_f);

            % integrated normal flux
            normal_flux = nout.' * we;

            wI_local(lf) = sum(normal_flux .* qw_f);

        end

        % MFD discrete inner product
        dw = wI_local - wh_local;
        AE = A_mfd_local{e};

        eW2 = eW2 + dw.' * AE * dw;
        nW2 = nW2 + wI_local.' * AE * wI_local;

    end

    % RELATIVE CELL ERRORS

    [err.displacement, err.displacementIsRelative] = relOrAbs(eU2, nU2, nU2);

    [err.divergence, err.divergenceIsRelative] = relOrAbs(eDiv2, nDiv2, nPi2);

    [err.projection, err.projectionIsRelative] = relOrAbs(ePi2, nPi2, nPi2);

    [err.pressure, err.pressureIsRelative] = relOrAbs(eP2, nP2, nP2);

    [err.flux, err.fluxIsRelative] = relOrAbs(eW2, nW2, nW2);

    % 6. TRACTION-BASED STRESS ERROR
    % E_sig^2 = sum_f h_f kappa int_f |(sigma - sigma_h)n|^2 df

    eS2 = 0;
    nS2 = 0;

    for f = 1:nFaces

        qp = face_struct(f).quad_points;
        qw = face_struct(f).quad_weights(:).';

        xf = face_struct(f).center(:);
        Af = face_struct(f).area;

        g = face_global_geom(f);

        % Face diameter
        P = V3(face_struct(f).verts(:), :);

        sq = sum(P.^2, 2);

        hf = sqrt(max(max(sq + sq.' - 2*(P*P.'), 0), [], 'all'));

        % kappa = 1/2 tr(C^{-1})
        e1 = face_struct(f).cells(1);
        kappa = 0.5 * trace(cell_struct(e1).Cinv);

        % Discrete traction
        c = sigma_h(6*(f-1) + (1:6));

        Th = zeros(3, numel(qw));

        for q = 1:numel(qw)

            Phi = evaluateFaceBasis(qp(:,q), xf, g.Qf, g.t1, g.t2, g.n, Af, g.m20, g.m02);
            Th(:,q) = Phi * c;

        end

        % Exact traction
        Se = problem.sigma(qp);

        Te = squeeze(sum(Se .* reshape(g.n, 1, 3, []), 2));

        if size(Te,2) ~= numel(qw)
            Te = reshape(Te, 3, []);
        end

        eS2 = eS2 + hf * kappa * sum(sum((Te - Th).^2, 1) .* qw);

        nS2 = nS2 + hf * kappa * sum(sum(Te.^2, 1) .* qw);

    end

    [err.stress, err.stressIsRelative] = relOrAbs(eS2, nS2, nS2);

end

% HELPER FUNCTIONS

% RELATIVE / ABSOLUTE ERROR
function [e, isRel] = relOrAbs(num, den, scale)

    if den > 1e-20 * max(scale, realmin)

        e = sqrt(num / den);
        isRel = true;

    else

        e = sqrt(num);
        isRel = false;

    end

end

% GLOBAL -> LOCAL STRESS DOFs
function sl = localStressDofs(e, cell_struct, face_struct, face_global_geom, B_local, sigma_h)

    fl = cell_struct(e).faces(:).';

    sl = zeros(6*numel(fl), 1);

    for lf = 1:numel(fl)

        f = fl(lf);

        Tf = computeFaceTransformation(f, B_local{e}.geom(lf), face_global_geom(f), face_struct);

        sl(6*(lf-1) + (1:6)) = diag(Tf) .* sigma_h(6*(f-1) + (1:6));

    end

end