% clear; clc; close all
addpath('geometry', 'operators', 'assembly', 'solvers', 'io', 'tests', genpath('FACTORIZE'), 'meshes');
% setupMRSTAuto();

%% STEP 1: LOAD MESH AND CHOOSE SYSTEM (mixed or hybrid), SOLVER TYPE (direct or iterative)
fprintf('\n[# MESH LOADED: unit cube model]\n');

meshName = 'unit_cube';

N = 8;
G = cartGrid([N N N],[1 1 1]);
G = computeGeometry(G);
[cell_struct, face_struct, V3, cells3D] = MRSTGridConvert(G);

%% STEP 2: SET UP GEOMETRY
% Lamé parameters (currently assuming homogeneous isotropic elasticity)
lambda = 1;
mu = 1;
alpha = 1;
s0 = 1;

cell_struct = assignCellProperties(cell_struct, V3, cells3D, lambda, mu, alpha, s0);
face_struct = assignFaceProperties(face_struct); 

cell_struct = computeVolumeQuadrature(cell_struct, face_struct, V3);
face_struct = computeFaceQuadrature(face_struct, V3);

nStressDofs = 6*numel(face_struct);
nDispDofs = 6*numel(cell_struct);

%% STEP 3: (LOCAL LEVEL) COMPUTE VEM MATRICES
fprintf('\n[# COMPUTE LOCAL OPERATORS]\n');

% INERTIA MATRIX (quadrature approach, Tonon approach)
A_local = computeAE_Tonon(cell_struct, face_struct, V3);

% DIVERGENCE MATRIX (quadrature approach, closed form approach from notes)
B_local = computeBE(cell_struct, face_struct, V3);
face_global_geom = computeGlobalFaceGeometry(face_struct, V3);

% PROJECTION MATRIX
D_local = computeDE(cell_struct, A_local, B_local);
Bproj = computeBproj(cell_struct, face_struct, D_local, B_local);
P_local = computeProjectionMatrix(cell_struct, Bproj);

% BIOT STRESS-PRESSURE COUPLING
D_biot_local = computeBiotCoupling(cell_struct, P_local);

K_cons = computeConsistency(cell_struct, P_local);
K_stab = computeStabilization(cell_struct, face_struct, B_local, P_local);

% FLOW PART (MIXED MFD)
A_mfd_local = computeMFDInnerProduct(cell_struct, face_struct);
E_mfd_local = computeMFDCoupling(cell_struct);
M_biot_local = computeBiotStorage(cell_struct);

%% STEP 4: GLOBAL ASSEMBLY
fprintf('\n[# ASSEMBLE GLOBAL MATRICES]\n');

[K_global, B_global, D_global, A_global, E_global, M_global] = assembleBiotGlobalMatrices(cell_struct, face_struct, face_global_geom, B_local, K_cons, K_stab, D_biot_local, A_mfd_local, E_mfd_local, M_biot_local);


%% STEP 5: TEMPORAL CONVERGENCE
T = 1.0;

% reference time step
dt_ref = 1/512;

% time steps for convergence test
dt_list = [1/2, 1/4, 1/8, 1/16];

%% 5.1 REFERENCE SOLUTION
fprintf('\n[# REFERENCE SOLUTION: dt = %.6f]\n', dt_ref);

[sigma_ref, u_ref, w_ref, p_ref] = solveBiotTime(T, dt_ref, cell_struct, face_struct, face_global_geom, K_global, B_global, D_global, A_global, E_global, M_global);


%% 5.2 TEMPORAL CONVERGENCE RUNS
nDt = numel(dt_list);

err_u = zeros(nDt,1);
% err_div = zeros(nDt,1);
err_proj = zeros(nDt,1);
err_stress = zeros(nDt,1);
err_p = zeros(nDt,1);
err_w = zeros(nDt,1);

for k = 1:nDt

    dt = dt_list(k);

    fprintf('\n[# TEMPORAL RUN: dt = %.6f]\n', dt);

    [sigma_h, u_h, w_h, p_h] = solveBiotTime(T, dt, cell_struct, face_struct, face_global_geom, K_global, B_global, D_global, A_global, E_global, M_global);

    % compare with reference solution at T
    err = computeBiotTemporalErrors(cell_struct, face_struct, face_global_geom, V3, B_local, D_local, P_local, A_mfd_local, sigma_h, u_h, w_h, p_h, sigma_ref, u_ref, w_ref, p_ref);

    err_u(k) = err.displacement;
    % err_div(k) = err.divergence;
    err_proj(k) = err.projection;
    err_stress(k) = err.stress;
    err_p(k) = err.pressure;
    err_w(k) = err.flux;

end


%% STEP 6: TEMPORAL CONVERGENCE RATES
rate_u = log(err_u(1:end-1)./err_u(2:end)) / log(2);
% rate_div = log(err_div(1:end-1)./err_div(2:end)) / log(2);
rate_proj = log(err_proj(1:end-1)./err_proj(2:end)) / log(2);
rate_stress = log(err_stress(1:end-1)./err_stress(2:end)) / log(2);
rate_p = log(err_p(1:end-1)./err_p(2:end)) / log(2);
rate_w = log(err_w(1:end-1)./err_w(2:end)) / log(2);

fprintf('\n============================================================\n');
fprintf('TEMPORAL CONVERGENCE (N = %d, T = %.2f)\n', N, T);
fprintf('Reference dt = %.6f\n', dt_ref);
fprintf('============================================================\n');

fprintf('\nErrors:\n');
fprintf('dt          Eu           Eproj       Esig        Ep          Ew\n');

for k = 1:nDt
    fprintf('%.6f   %.3e   %.3e   %.3e   %.3e   %.3e\n', dt_list(k), err_u(k), err_proj(k), err_stress(k), err_p(k), err_w(k));
end

fprintf('\nRates:\n');
fprintf('dt ratio    Eu       Eproj    Esig     Ep       Ew\n');

for k = 1:nDt-1
    fprintf('%g -> %g   %.3f    %.3f    %.3f    %.3f    %.3f\n', dt_list(k), dt_list(k+1), rate_u(k), rate_proj(k), rate_stress(k), rate_p(k), rate_w(k));
end

fprintf('============================================================\n');


%% STEP 7: PLOT
figure;

loglog(1./dt_list, err_u,      '-o', 'LineWidth', 1.5, 'MarkerSize', 7); 
hold on;
loglog(1./dt_list, err_proj,   '-s', 'LineWidth', 1.5, 'MarkerSize', 7);
loglog(1./dt_list, err_stress, '-^', 'LineWidth', 1.5, 'MarkerSize', 7);
loglog(1./dt_list, err_p,      '-d', 'LineWidth', 1.5, 'MarkerSize', 7);
loglog(1./dt_list, err_w,      '-x', 'LineWidth', 1.5, 'MarkerSize', 8);

% O(dt) reference line
ref1 = err_p(1) * (dt_list ./ dt_list(1));
loglog(1./dt_list, ref1, '--', 'LineWidth', 1.5);

grid on;
box on;

xlabel('1/\Delta t');
ylabel('Relative temporal error');
title('Temporal Convergence: Mixed VEM - Mixed MFD Biot');

legend('$E_u$', '$E_{\sigma,\Pi}$', '$E_{\sigma}$', '$E_p$', '$E_w$', '$O(\Delta t)$', 'Interpreter','latex', 'Location','best');
set(gca,'FontSize',12);


%% LOCAL FUNCTION: BIOT TIME STEPPING
function [sigma_h, u_h, w_h, p_h] = solveBiotTime(T, dt, cell_struct, face_struct, face_global_geom, K_global, B_global, D_global, A_global, E_global, M_global)

    nSteps = round(T/dt);

    nStressDofs = size(K_global,1);
    nDispDofs = size(B_global,1);
    nFluxDofs = size(A_global,1);
    nPressureDofs = size(M_global,1);

    % backward Euler system matrix
    LHS = [K_global, B_global', sparse(nStressDofs,nFluxDofs), D_global';
           -B_global, sparse(nDispDofs,nDispDofs), sparse(nDispDofs,nFluxDofs), sparse(nDispDofs,nPressureDofs);
           sparse(nFluxDofs,nStressDofs), sparse(nFluxDofs,nDispDofs), A_global, E_global';
           D_global, sparse(nPressureDofs,nDispDofs), -dt*E_global, M_global ];

    % exact discrete initial condition at t = 0
    [sigma_old, p_old] = initializeBiot(cell_struct, face_struct, face_global_geom);

    % time stepping
    for n = 1:nSteps

        t = n * dt;

        [f_sigma, f_u, f_w, f_p] = assembleBiotRHS(cell_struct, face_struct, face_global_geom, t);

        RHS = [ ...
            f_sigma;
            f_u;
            f_w;
            D_global*sigma_old ...
            + M_global*p_old ...
            + dt*f_p ];

        solution = LHS \ RHS;

        i1 = nStressDofs;
        i2 = i1 + nDispDofs;
        i3 = i2 + nFluxDofs;

        sigma_h = solution(1:i1);
        u_h = solution(i1+1:i2);
        w_h = solution(i2+1:i3);
        p_h = solution(i3+1:end);

        sigma_old = sigma_h;
        p_old = p_h;

    end

end

function err = computeBiotTemporalErrors(cell_struct, face_struct, face_global_geom, V3, B_local, D_local, P_local, A_mfd_local, sigma_h, u_h, w_h, p_h, sigma_ref, u_ref, w_ref, p_ref)

    nCells = numel(cell_struct);
    nFaces = numel(face_struct);

    eU2 = 0;   nU2 = 0;
    eDiv2 = 0; nDiv2 = 0;
    ePi2 = 0;  nPi2 = 0;
    eP2 = 0;   nP2 = 0;
    eW2 = 0;   nW2 = 0;

    s2 = 1/sqrt(2);

    %% cell-based errors
    for e = 1:nCells

        qp = cell_struct(e).quad_points;
        qw = cell_struct(e).quad_weights(:).';

        xE = cell_struct(e).center(:);
        r = qp - xE;

        % local stress DOFs
        fl = cell_struct(e).faces(:).';
        nf = numel(fl);

        sl_h = zeros(6*nf,1);
        sl_ref = zeros(6*nf,1);

        for lf = 1:nf

            f = fl(lf);

            Tf = computeFaceTransformation(f, B_local{e}.geom(lf), face_global_geom(f), face_struct);

            ids = 6*(f-1) + (1:6);
            lid = 6*(lf-1) + (1:6);

            sl_h(lid) = diag(Tf) .* sigma_h(ids);
            sl_ref(lid) = diag(Tf) .* sigma_ref(ids);

        end

        %% displacement
        idsU = 6 * (e-1) + (1:6);

        uh_coeff = u_h(idsU);
        uref_coeff = u_ref(idsU);

        uh = uh_coeff(1:3) + cross(repmat(uh_coeff(4:6),1,size(r,2)), r, 1);

        uref = uref_coeff(1:3) + cross(repmat(uref_coeff(4:6),1,size(r,2)), r, 1);

        eU2 = eU2 + sum(sum((uh-uref).^2,1).*qw);
        nU2 = nU2 + sum(sum(uref.^2,1).*qw);


        %% stress divergence
        dh = D_local{e} * sl_h;
        dref = D_local{e} * sl_ref;

        divh = dh(1:3) + cross(repmat(dh(4:6),1,size(r,2)), r, 1);

        divref = dref(1:3) + cross(repmat(dref(4:6),1,size(r,2)), r, 1);

        eDiv2 = eDiv2 + sum(sum((divh-divref).^2,1).*qw);
        nDiv2 = nDiv2 + sum(sum(divref.^2,1).*qw);


        %% projected stress
        ch = P_local{e} * sl_h;
        cref = P_local{e} * sl_ref;

        Sh = [ ...
            ch(1),    s2*ch(6), s2*ch(5);
            s2*ch(6), ch(2),    s2*ch(4);
            s2*ch(5), s2*ch(4), ch(3)];

        Sref = [ ...
            cref(1),    s2*cref(6), s2*cref(5);
            s2*cref(6), cref(2),    s2*cref(4);
            s2*cref(5), s2*cref(4), cref(3)];

        volume = sum(qw);

        ePi2 = ePi2 + volume * sum((Sh-Sref).^2,'all');
        nPi2 = nPi2 + volume * sum(Sref.^2,'all');


        %% pressure
        dp = p_h(e) - p_ref(e);

        eP2 = eP2 + volume * dp^2;
        nP2 = nP2 + volume * p_ref(e)^2;


        %% flux
        face_ids = cell_struct(e).faces(:);
        orient = cell_struct(e).faces_orientation(:);

        wh   = orient .* w_h(face_ids);
        wref = orient .* w_ref(face_ids);

        dw = wh - wref;

        AE = A_mfd_local{e};

        eW2 = eW2 + dw' * AE * dw;
        nW2 = nW2 + wref' * AE * wref;

    end


    %% traction stress error
    eS2 = 0;
    nS2 = 0;

    for f = 1:nFaces

        qp = face_struct(f).quad_points;
        qw = face_struct(f).quad_weights(:).';

        xf = face_struct(f).center(:);
        Af = face_struct(f).area;

        g = face_global_geom(f);

        % face diameter
        P = V3(face_struct(f).verts(:),:);
        d2 = sum(P.^2,2) + sum(P.^2,2).' - 2*(P*P.');
        hf = sqrt(max(d2,[],'all'));

        e1 = face_struct(f).cells(1);
        kappa = 0.5 * trace(cell_struct(e1).Cinv);

        ids = 6 * (f-1) + (1:6);

        ch = sigma_h(ids);
        cref = sigma_ref(ids);

        for q = 1:numel(qw)

            Phi = evaluateFaceBasis(qp(:,q), xf, g.Qf, g.t1, g.t2, g.n, Af, g.m20, g.m02);

            th = Phi * ch;
            tref = Phi * cref;

            eS2 = eS2 + hf*kappa*qw(q)*sum((th-tref).^2);
            nS2 = nS2 + hf*kappa*qw(q)*sum(tref.^2);

        end

    end


    %% relative errors
    err.displacement = sqrt(eU2/nU2);
    % err.divergence = sqrt(eDiv2/nDiv2);
    err.projection = sqrt(ePi2/nPi2);
    err.stress = sqrt(eS2/nS2);
    err.pressure = sqrt(eP2/nP2);
    err.flux = sqrt(eW2/nW2);

end


