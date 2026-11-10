% clear; clc; close all
addpath('geometry', 'operators', 'assembly', 'solvers', 'io', 'tests', genpath('FACTORIZE'), 'meshes');
% setupMRSTAuto();

%% STEP 1: LOAD MESH
fprintf('\n[# MESH LOADED: unit cube model]\n');

val = 32;
Nx = val;
Ny = val;
Nz = val;

G = cartGrid([Nx Ny Nz],[1 1 1]);
G = computeGeometry(G);

[cell_struct, face_struct, V3, cells3D] = MRSTGridConvert(G);

%% STEP 2: SET UP GEOMETRY
% Lamé parameters (currently assuming homogeneous isotropic elasticity)
lambda = 1;
mu = 1;
alpha = 1;
s0 = 1;
xFrac = 0.5;
yRange = [0, 1.0];
zRange = [0, 1.0];

cell_struct = assignCellProperties(cell_struct, V3, cells3D, lambda, mu, alpha, s0);
face_struct = assignFaceProperties(face_struct); 

% assign DOFs over the fracture
[cell_struct, face_struct, fracture] = assignFractureProperties(cell_struct, face_struct, xFrac, yRange, zRange);

cell_struct = computeVolumeQuadrature(cell_struct, face_struct, V3);
face_struct = computeFaceQuadrature(face_struct, V3);

% DOFs assignment
% [cell_struct, face_struct, nStressDofs] = assignFractureStressDofs(cell_struct, face_struct, fracture);
fracture = buildFractureSurfaceMesh(fracture, face_struct, V3);
[cell_struct, face_struct, fracture, flowDofs] = assignFractureFlowDofs(cell_struct, face_struct, fracture);

%% STEP 3: COMPUTE LOCAL MFD OPERATORS
fprintf('\n[# COMPUTE LOCAL MFD OPERATORS]\n');

% Bulk Darcy flow
A_mfd_local = computeMFDInnerProduct(cell_struct, face_struct);
E_mfd_local = computeMFDCoupling(cell_struct);

% Fracture Darcy flow
AF_mfd_local = computeFractureMFDInnerProduct(fracture);
EF_mfd_local = computeFractureMFDCoupling(fracture);

%% STEP 4: COMPUTE INTERFACE OPERATORS
fprintf('\n[# COMPUTE INTERFACE OPERATORS]\n');

Tminus = 1;
Tplus = 1;

R = computeFractureRobinOperator(fracture, face_struct, flowDofs, Tminus, Tplus);
TF = computeFractureTransferOperator(fracture, cell_struct, face_struct, flowDofs);

%% STEP 5: GLOBAL FLOW ASSEMBLY
[A, flowBlocks] = assembleMixedDimensionalFlowSystem(cell_struct, fracture, flowDofs, A_mfd_local, E_mfd_local, AF_mfd_local, EF_mfd_local, R, TF);

% b = assembleFlowRHS(cell_struct, fracture, flowDofs);
exactFun = @(x) fullFractureExactData(x, xFrac);

b = assembleFlowRHS(cell_struct, fracture, flowDofs, face_struct, exactFun);

[Ared, bred, bc] = applyFlowBC(A, b, cell_struct, face_struct, fracture, flowDofs, exactFun);


%% STEP 6: SOLVE GLOBAL FLOW SYSTEM
fprintf('\n[# SOLVE MIXED-DIMENSIONAL FLOW SYSTEM]\n');

xred = Ared \ bred;

% Recover full solution including prescribed boundary DOFs
x = zeros(size(A, 1), 1);

x(bc.freeDofs) = xred;
x(bc.fixedDofs) = bc.fixedValues;

% Extract solution components
w = x(1:flowDofs.nW);

p = x(flowDofs.offset.p + (1:flowDofs.nP));

wF = x(flowDofs.offset.wF + (1:flowDofs.nWF));

pF = x(flowDofs.offset.pF + (1:flowDofs.nPF));


%% STEP 7: COMPARE WITH MANUFACTURED EXACT SOLUTION
fprintf('\n[# MANUFACTURED SOLUTION ERROR]\n');

% BULK PRESSURE ERROR
pExact = zeros(flowDofs.nP, 1);
bulkVolumes = zeros(flowDofs.nP, 1);

for e = 1:numel(cell_struct)

    xc = cell_struct(e).center(:);

    exact = fullFractureExactData(xc, xFrac);

    pExact(e) = exact.pBulk;
    bulkVolumes(e) = cell_struct(e).volume;
end

% Volume-weighted relative discrete L2 error
bulkError = sqrt(sum(bulkVolumes .* (p - pExact).^2));
bulkNorm = sqrt(sum(bulkVolumes .* pExact.^2));

relBulkPressureError = bulkError / max(bulkNorm, eps);

% FRACTURE PRESSURE ERROR
pFExact = zeros(flowDofs.nPF, 1);
fractureAreas = zeros(flowDofs.nPF, 1);

for k = 1:fracture.nFaces

    f = fracture.faces(k);
    dof = fracture.pressure_dofs(k);

    xf = face_struct(f).center(:);

    exact = fullFractureExactData(xf, xFrac);

    pFExact(dof) = exact.pFrac;
    fractureAreas(dof) = fracture.face_areas(k);
end

% Area-weighted relative discrete L2 error
fractureError = sqrt(sum(fractureAreas .* (pF - pFExact).^2));
fractureNorm = sqrt(sum(fractureAreas .* pFExact.^2));

relFracturePressureError = fractureError / max(fractureNorm, eps);

fprintf('  Bulk pressure absolute error: %.6e\n', bulkError);
fprintf('  Bulk pressure relative error: %.6e\n', relBulkPressureError);

fprintf('  Fracture pressure absolute error: %.6e\n', fractureError);
fprintf('  Fracture pressure relative error: %.6e\n', relFracturePressureError);

%% STEP 8: FLUX CONVERGENCE ERROR
fprintf('\n[# MANUFACTURED FLUX ERROR]\n');

% Exact bulk flux DOFs
wExact = zeros(flowDofs.nW, 1);

for f = 1:numel(face_struct)

    xf = face_struct(f).center(:);
    Af = face_struct(f).area;

    n = face_struct(f).normal(:);
    n = n / norm(n);

    exact = exactFun(xf);

    if face_struct(f).is_fracture

        % Original DOF: minus-side cell
        dofMinus = face_struct(f).flow_dof;

        % Duplicated DOF: plus-side cell
        dofPlus = face_struct(f).flow_dof_plus;

        wExact(dofMinus) = Af * dot(exact.wBulkMinus, n);

        wExact(dofPlus) = Af * dot(exact.wBulkPlus, n);

    else

        % Ordinary bulk face
        wExact(face_struct(f).flow_dof) = Af * dot(exact.wBulk, n);

    end
end

% Exact fracture edge flux DOFs
wFExact = zeros(flowDofs.nWF, 1);

for e = 1:fracture.nEdges

    xe = fracture.edge_centers(e,:)';
    Le = fracture.edge_lengths(e);

    exact = exactFun(xe);

    nu = fracture.edge_normals(e,:)';
    nu = nu / norm(nu);

    wFExact(fracture.edge_flux_dofs(e)) = Le * dot(exact.wFrac, nu);
end

% Extract pure bulk/fracture MFD inner products
W = flowBlocks.W;
WF = flowBlocks.WF;

% MFD norm errors
dw = w - wExact;
dwF = wF - wFExact;

bulkFluxError = sqrt(max(dw' * W * dw, 0));
fractureFluxError = sqrt(max(dwF' * WF * dwF, 0));

bulkFluxNorm = sqrt(max(wExact' * W * wExact, 0));
fractureFluxNorm = sqrt(max(wFExact' * WF * wFExact, 0));

relBulkFluxError = bulkFluxError / max(bulkFluxNorm, eps);
relFractureFluxError = fractureFluxError / max(fractureFluxNorm, eps);

fprintf('  Bulk flux absolute error: %.6e\n', bulkFluxError);
fprintf('  Bulk flux relative error: %.6e\n', relBulkFluxError);

fprintf('  Fracture flux absolute error: %.6e\n', fractureFluxError);
fprintf('  Fracture flux relative error: %.6e\n', relFractureFluxError);





%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function exact = fullFractureExactData(x, xFrac)
% exact manufactured solution for a full planar fracture.

    x = x(:);

    assert(numel(x) == 3, 'Expected a 3D coordinate.');
    assert(all(isfinite(x)), 'Invalid coordinate.');

    y = x(2);
    z = x(3);

    d = abs(x(1) - xFrac);

    a = y * (1 - y);
    c = z * (1 - z);

    omega = a^2 * c^2;

    lapOmegaY = 2 * ((1 - 2*y)^2 - 2*a) * c^2;
    lapOmegaZ = 2 * ((1 - 2*z)^2 - 2*c) * a^2;

    lapOmega = lapOmegaY + lapOmegaZ;

    exact.pBulk = d^(5/2) + omega * d;
    exact.pFrac = -omega;

    exact.lambdaMinus = omega;
    exact.lambdaPlus = omega;

    exact.fBulk = -(15/4) * sqrt(d) - d * lapOmega;
    exact.fFrac = lapOmega - 2 * omega;

    exact.omega = omega;
    exact.d = d;

    %% EXACT BULK DARCY FLUX

    domegaDy = 2 * a * (1 - 2*y) * c^2;
    domegaDz = 2 * c * (1 - 2*z) * a^2;

    wx = -sign(x(1) - xFrac) * ((5/2) * d^(3/2) + omega);

    wy = -d * domegaDy;
    wz = -d * domegaDz;

    exact.wBulk = [wx; wy; wz];

    % One-sided traces on the fracture
    exact.wBulkMinus = [ omega; -d*domegaDy; -d*domegaDz ];
    exact.wBulkPlus  = [-omega; -d*domegaDy; -d*domegaDz ];

    %% EXACT FRACTURE DARCY FLUX
    exact.wFrac = [0; domegaDy; domegaDz];
end
