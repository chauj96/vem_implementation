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


%% STEP 5: SOLVE THE BIOT PROBLEM
tau = 4/3;

nSteps = 1600;
dt = tau/nSteps;
T = tau;

fprintf('\n[# CONSOLIDATION SIMULATION]\n');
fprintf('Final time: %.4f, dt: %.4f\n', T, dt);

[sigma_h, u_h, w_h, p_h, p_history, t_history] = solveBiotTime(T, dt, cell_struct, face_struct, face_global_geom, K_global, B_global, D_global, A_global, E_global, M_global);

fprintf('Pressure history size: %d x %d\n', size(p_history,1), size(p_history,2));
fprintf('Stored time range: [%.4f, %.4f]\n', t_history(1), t_history(end));
fprintf('Final pressure difference: %.3e\n', norm(p_history(:,end) - p_h));

%% STEP 6: COMPARE WITH ANALYTICAL SOLUTION

% Parameters
H = lambda + 2*mu;
varpi = 1;
Kzz = 1;

p_init = alpha*varpi / (s0*H + alpha^2);
cf = Kzz / (s0 + alpha^2/H);

% Normalized time
tau = 1/cf;
tbar_targets = [0.01, 0.05, 0.1, 0.2, 0.5, 1];

% Normalized depth
z = arrayfun(@(c) c.center(3), cell_struct).';

[z_levels, ~, depthIdx] = unique(z);
z_fine = linspace(0, 1, 300)';

% Plot
colors = lines(numel(tbar_targets));

figure;
hold on;

for k = 1:numel(tbar_targets)

    tbar = tbar_targets(k);

    % NUMERICAL SOLUTION
    t_target = tbar * tau;
    [~, timeIdx] = min(abs(t_history - t_target));

    p_num = p_history(:,timeIdx) / p_init;
    p_profile = accumarray(depthIdx, p_num, [], @mean);

    % ANALYTICAL SOLUTION
    p_exact = zeros(size(z_fine));

    for n = 0:100
        kn = (2*n+1)*pi/2;
        p_exact = p_exact + (4/((2*n+1)*pi)) * sin(kn*z_fine) * exp(-kn^2*tbar);
    end

    plot(p_exact, z_fine, '-', 'Color', colors(k,:), 'LineWidth', 1.8, 'DisplayName', sprintf('Analytical, tbar = %.3g', tbar));
    plot(p_profile, z_levels, 'o', 'Color', colors(k,:), 'MarkerSize', 6, 'LineWidth', 1.2, 'DisplayName', sprintf('Mixed VEM-MFD, tbar = %.3g', tbar));

end

% Figure settings
set(gca, 'YDir', 'reverse');

xlabel('Normalized pressure, p / p_{init}');
ylabel('Normalized depth, z / h');
title('Consolidation of a Soil Layer');

xlim([0 1.05]);
ylim([0 1]);

grid on;
legend('Location', 'eastoutside');


%% LOCAL FUNCTION: BIOT TIME STEPPING
function [sigma_h, u_h, w_h, p_h, p_history, t_history] = solveBiotTime(T, dt, cell_struct, face_struct, face_global_geom, K_global, B_global, D_global, A_global, E_global, M_global)

    nSteps = round(T/dt);

    nStressDofs = size(K_global,1);
    nDispDofs = size(B_global,1);
    nFluxDofs = size(A_global,1);
    nPressureDofs = size(M_global,1);

    % global DOF offsets
    fluxOffset = nStressDofs + nDispDofs;

    % backward Euler system matrix
    LHS = [ ...
        K_global, B_global', sparse(nStressDofs,nFluxDofs), D_global';
        -B_global, sparse(nDispDofs,nDispDofs), sparse(nDispDofs,nFluxDofs), sparse(nDispDofs,nPressureDofs);
        sparse(nFluxDofs,nStressDofs), sparse(nFluxDofs,nDispDofs), A_global, E_global';
        D_global, sparse(nPressureDofs,nDispDofs), -dt*E_global, M_global ];

    nTotalDofs = size(LHS,1);

    %% BOUNDARY CONDITIONS

    fixedDofs = [];
    fixedValues = [];

    varpi = 1;
    ez = [0; 0; 1];

    for f = 1:numel(face_struct)

        face = face_struct(f);
        geom = face_global_geom(f);

        stressIds = 6*(f-1) + (1:6);
        fluxId = fluxOffset + f;

        if face.is_top

            % prescribed total traction: sigma*n_out = varpi*ez
            e = face.cells(1);

            outwardDirection = face.center(:) - cell_struct(e).center(:);
            normalSign = sign(dot(geom.n, outwardDirection));

            tractionGlobal = normalSign * varpi * ez;

            Af = face.area;

            tractionCoeffs = zeros(6,1);

            tractionCoeffs(1) = Af * dot(geom.t1,tractionGlobal);
            tractionCoeffs(2) = Af * dot(geom.t2, tractionGlobal);
            tractionCoeffs(3) = Af * dot(geom.n, tractionGlobal);

            % modes 4, 5, 6 vanish for constant traction

            fixedDofs = [fixedDofs; stressIds(:)];
            fixedValues = [fixedValues; tractionCoeffs];

        elseif face.is_bottom || face.is_side

            % zero tangential traction:
            % modes 1, 2, 4
            tangentialIds = stressIds([1,2,4]);

            fixedDofs = [fixedDofs; tangentialIds(:)];
            fixedValues = [fixedValues; zeros(3,1)];

            % impermeable boundary: w*n = 0
            fixedDofs = [fixedDofs; fluxId];
            fixedValues = [fixedValues; 0];

        end

    end

    %% REDUCED SYSTEM

    allDofs = (1:nTotalDofs).';

    freeDofs = setdiff(allDofs, fixedDofs);

    LHS_free = LHS(freeDofs, freeDofs);
    LHS_fixed = LHS(freeDofs, fixedDofs);

    % factorize once because LHS is constant in time
    solver = decomposition(LHS_free, 'lu');

    %% INITIAL CONDITION

    [sigma_old, p_old] = initializeBiot(cell_struct, face_struct, face_global_geom);

    p_history = zeros(nPressureDofs, nSteps+1);
    t_history = (0:nSteps)' * dt;
    
    % initial pressure at t = 0
    p_history(:,1) = p_old;

    %% TIME STEPPING

    for n = 1:nSteps

        t = n*dt;

        [f_sigma, f_u, f_w, f_p] = assembleBiotRHS(cell_struct, face_struct, face_global_geom, t);

        RHS = [ ...
            f_sigma;
            f_u;
            f_w;
            D_global*sigma_old + M_global*p_old + dt*f_p ];

        % account for prescribed boundary values
        RHS_free = RHS(freeDofs) - LHS_fixed*fixedValues;

        % solve reduced system
        solution = zeros(nTotalDofs,1);

        solution(fixedDofs) = fixedValues;
        solution(freeDofs) = solver \ RHS_free;

        %% EXTRACT SOLUTION

        i1 = nStressDofs;
        i2 = i1 + nDispDofs;
        i3 = i2 + nFluxDofs;

        sigma_h = solution(1:i1);
        u_h = solution(i1+1:i2);
        w_h = solution(i2+1:i3);
        p_h = solution(i3+1:end);

        p_history(:,n+1) = p_h;
        
        sigma_old = sigma_h;
        p_old = p_h;

    end

end