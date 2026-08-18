function err = computeErrorsDassi(cell_struct, face_struct, face_global_geom, V3, ...
                                 B_local, D_local, P_local, sigma_h, u_h, problem)
% error norms of Dassi, Lovadina & Visinoni (2020), section 5.1
%
%   E_u      = || u - u_h ||_0                                    (full L2, not a
%                                                                  point value)
%   E_sigdiv = || div sigma - div sigma_h ||_0
%   E_sigPi  = || sigma - Pi_E sigma_h ||_0
%   E_sig    = ( sum_f h_f int_f kappa |(sigma - sigma_h) n|^2 df )^1/2
%
% note on E_sig: the discrete traction on a face is the FULL element of T_h(f),
% i.e. all six basis functions including the two normal-moment and one in-plane
% rotation terms. it is weighted by the face diameter h_f, not the face area,
% and by kappa = 1/2 tr(C^-1). the paper expects these to converge at rate ~ h
%
% u_h|_E = alpha + omega x (x - x_E) is evaluated as a field, and
% div sigma_h|_E = alpha_E + omega_E x (x - x_E) comes from D_local.
%
% problem.f holds div(sigma) in this repository's sign convention, i.e. the
% discrete equation is B sigma = int f . v. the paper writes -div sigma = f, so
% its loading term is the negative of this one

    nCells = numel(cell_struct);
    nFaces = numel(face_struct);

    eU2 = 0; eDiv2 = 0; ePi2 = 0;
    nU2 = 0; nDiv2 = 0; nPi2 = 0;

    s = 1/sqrt(2);

    for e = 1:nCells

        qp = cell_struct(e).quad_points;
        qw = cell_struct(e).quad_weights(:).';
        xE = cell_struct(e).center(:);
        r = qp - xE;

        % local stress DOFs in the element's own face convention
        sl = localStressDofs(e, cell_struct, face_struct, face_global_geom, B_local, sigma_h);

        %% E_u : u_h = alpha + omega x r
        a = u_h(6*(e-1) + (1:3));
        w = u_h(6*(e-1) + (4:6));

        uh = a + cross(repmat(w,1,size(r,2)), r, 1);
        ue = problem.u(qp);

        eU2 = eU2 + sum(sum((ue - uh).^2, 1) .* qw);
        nU2 = nU2 + sum(sum(ue.^2, 1) .* qw);

        %% E_sigdiv : div sigma_h = alphaE + omegaE x r
        d = D_local{e} * sl;
        dh = d(1:3) + cross(repmat(d(4:6),1,size(r,2)), r, 1);
        de = problem.f(qp);              % = div(sigma), see the header note

        eDiv2 = eDiv2 + sum(sum((de - dh).^2, 1) .* qw);
        nDiv2 = nDiv2 + sum(sum(de.^2, 1) .* qw);

        %% E_sigPi : Pi_E sigma_h, a constant symmetric tensor
        c = P_local{e} * sl;
        Sh = [   c(1), s*c(6), s*c(5); ...
              s*c(6),    c(2), s*c(4); ...
              s*c(5), s*c(4),    c(3)];

        Se = problem.sigma(qp);
        dS = Se - Sh;

        ePi2 = ePi2 + sum(reshape(sum(sum(dS.^2,1),2), 1, []) .* qw);
        nPi2 = nPi2 + sum(reshape(sum(sum(Se.^2,1),2), 1, []) .* qw);

    end

    % relative where the reference is meaningful, absolute where it vanishes.
    % div(sigma) is identically zero for an unloaded body, so normalising by it
    % would report a meaningless blow-up rather than an error
    [err.displacement, err.displacementIsRelative] = relOrAbs(eU2, nU2, nU2);
    [err.divergence,   err.divergenceIsRelative]   = relOrAbs(eDiv2, nDiv2, nPi2);
    [err.projection,   err.projectionIsRelative]   = relOrAbs(ePi2, nPi2, nPi2);

    %% E_sig : full traction error on each face, weighted by h_f and kappa

    eS2 = 0; nS2 = 0;

    for f = 1:nFaces

        qp = face_struct(f).quad_points;
        qw = face_struct(f).quad_weights(:).';

        xf = face_struct(f).center(:);
        Af = face_struct(f).area;
        g = face_global_geom(f);

        % face diameter
        P = V3(face_struct(f).verts(:), :);
        sq = sum(P.^2, 2);
        hf = sqrt(max(max(sq + sq.' - 2*(P*P.'), 0), [], 'all'));

        % kappa = 1/2 tr(C^-1) of an adjacent cell
        e1 = face_struct(f).cells(1);
        kappa = 0.5 * trace(cell_struct(e1).Cinv);

        % discrete traction: the FULL element of T_h(f), all six basis functions
        c = sigma_h(6*(f-1) + (1:6));

        Th = zeros(3, numel(qw));
        for q = 1:numel(qw)
            Phi = evaluateFaceBasis(qp(:,q), xf, g.Qf, g.t1, g.t2, g.n, Af, g.m20, g.m02);
            Th(:,q) = Phi * c;
        end

        Se = problem.sigma(qp);
        Te = squeeze(sum(Se .* reshape(g.n, 1, 3, []), 2));
        if size(Te,2) ~= numel(qw), Te = reshape(Te, 3, []); end

        eS2 = eS2 + hf * kappa * sum(sum((Te - Th).^2, 1) .* qw);
        nS2 = nS2 + hf * kappa * sum(sum(Te.^2, 1) .* qw);

    end

    [err.stress, err.stressIsRelative] = relOrAbs(eS2, nS2, nS2);

end

function [e, isRel] = relOrAbs(num, den, scale)
% relative error when the reference norm is significant against scale,
% absolute error otherwise

    if den > 1e-20 * max(scale, realmin)
        e = sqrt(num / den);
        isRel = true;
    else
        e = sqrt(num);
        isRel = false;
    end

end

function sl = localStressDofs(e, cell_struct, face_struct, face_global_geom, B_local, sigma_h)
% sigma_h restricted to element e, converted to the element's own face
% convention: x_local = T_E x_global

    fl = cell_struct(e).faces(:).';
    sl = zeros(6*numel(fl), 1);

    for lf = 1:numel(fl)
        f = fl(lf);
        Tf = computeFaceTransformation(f, B_local{e}.geom(lf), face_global_geom(f), face_struct);
        sl(6*(lf-1) + (1:6)) = diag(Tf) .* sigma_h(6*(f-1) + (1:6));
    end

end
