% clear; clc; close all
addpath('geometry', 'operators', 'assembly', 'solvers', 'io', 'tests', genpath('FACTORIZE'), 'meshes');
% setupMRSTAuto();

%% STEP 1: LOAD MESH
fprintf('\n[# MESH LOADED: unit cube model]\n');

Lx = 1;
Ly = 1;
Lz = 1;

val = 32;

Nx = val;
Ny = val;
Nz = val;
meshType = 'tetra';

xFrac = 0.5;
yRange = [0, 1];
zRange = [0, 1];

[G, fractureInfo] = unitCubeFracture(Nx, Ny, Nz, xFrac, yRange, zRange, meshType);
[cell_struct, face_struct, V3, cells3D] = MRSTGridConvert(G);

% meshName = sprintf('unit_cube_%s_%d', meshType, val);
% save([meshName, '.mat'], 'cell_struct', 'face_struct', 'V3', 'cells3D', 'fractureInfo', 'Lx', 'Ly', 'Lz', '-v7');

figure('Color','w');
plotGrid(G, 'FaceColor','none', 'EdgeColor',[0.5,0.5,0.5]);
axis equal tight;
view(35,25);
xlabel('x');
ylabel('y');
zlabel('z');
title(sprintf('%s mesh: %d cells', meshType, G.cells.num));

%% STEP 2: SET UP GEOMETRY
% Lamé parameters (currently assuming homogeneous isotropic elasticity)
lambda = 1;
mu = 1;
alpha = 1;
s0 = 1;

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

% bulk pressure errors
bulkErrorCentroid = 0;
bulkErrorProjected = 0;
bulkErrorL2 = 0;

for e = 1:numel(cell_struct)

    xc = cell_struct(e).center(:);
    Ve = cell_struct(e).volume;

    qp = cell_struct(e).quad_points;
    qw = cell_struct(e).quad_weights;

    pe = p(e);

    % centroid-based error
    exact = exactFun(xc);

    bulkErrorCentroid = bulkErrorCentroid + Ve*(pe-exact.pBulk)^2;

    % quadrature-based integration
    pIntegral = 0;

    for q = 1:numel(qw)

        exact = exactFun(qp(:,q));

        pIntegral = pIntegral + qw(q)*exact.pBulk;

        bulkErrorL2 = bulkErrorL2 + qw(q)*(pe-exact.pBulk)^2;

    end

    % cell-average projection
    pAverage = pIntegral/Ve;

    bulkErrorProjected = bulkErrorProjected + Ve*(pe-pAverage)^2;

end

bulkErrorCentroid = sqrt(bulkErrorCentroid);
bulkErrorProjected = sqrt(bulkErrorProjected);
bulkErrorL2 = sqrt(bulkErrorL2);


% fracture pressure errors
fractureErrorCentroid = 0;
fractureErrorProjected = 0;
fractureErrorL2 = 0;

for k = 1:fracture.nFaces

    f = fracture.faces(k);
    dof = fracture.pressure_dofs(k);

    xf = face_struct(f).center(:);
    Af = fracture.face_areas(k);

    qp = face_struct(f).quad_points;
    qw = face_struct(f).quad_weights;

    pf = pF(dof);

    % centroid-based error
    exact = exactFun(xf);

    fractureErrorCentroid = fractureErrorCentroid + Af*(pf-exact.pFrac)^2;

    % quadrature-based integration
    pIntegral = 0;

    for q = 1:numel(qw)

        exact = exactFun(qp(:,q));

        pIntegral = pIntegral + qw(q)*exact.pFrac;

        fractureErrorL2 = fractureErrorL2 + qw(q)*(pf-exact.pFrac)^2;

    end

    % face-average projection
    pAverage = pIntegral/Af;

    fractureErrorProjected = fractureErrorProjected + Af*(pf-pAverage)^2;

end

fractureErrorCentroid = sqrt(fractureErrorCentroid);
fractureErrorProjected = sqrt(fractureErrorProjected);
fractureErrorL2 = sqrt(fractureErrorL2);


% print pressure errors
fprintf('\n[# PRESSURE ERRORS: CENTROID-BASED]\n');
fprintf('  Bulk pressure error:      %.6e\n', bulkErrorCentroid);
fprintf('  Fracture pressure error:  %.6e\n', fractureErrorCentroid);

fprintf('\n[# PRESSURE ERRORS: L2 PROJECTION]\n');
fprintf('  Bulk pressure error:      %.6e\n', bulkErrorProjected);
fprintf('  Fracture pressure error:  %.6e\n', fractureErrorProjected);

fprintf('\n[# PRESSURE ERRORS: CONTINUOUS L2]\n');
fprintf('  Bulk pressure error:      %.6e\n', bulkErrorL2);
fprintf('  Fracture pressure error:  %.6e\n', fractureErrorL2);

%% STEP 8: FLUX CONVERGENCE ERROR
fprintf('\n[# MANUFACTURED FLUX ERROR]\n');

% exact bulk flux DOFs
wExact = zeros(flowDofs.nW, 1);

for f = 1:numel(face_struct)

    xf = face_struct(f).center(:);
    Af = face_struct(f).area;

    n = face_struct(f).normal(:);
    n = n / norm(n);

    exact = exactFun(xf);

    if face_struct(f).is_fracture

        % minus-side cell
        dofMinus = face_struct(f).flow_dof;

        % plus-side cell
        dofPlus = face_struct(f).flow_dof_plus;

        wExact(dofMinus) = Af * dot(exact.wBulkMinus, n);

        wExact(dofPlus) = Af * dot(exact.wBulkPlus, n);

    else
        wExact(face_struct(f).flow_dof) = Af * dot(exact.wBulk, n);

    end
end

wFExact = zeros(flowDofs.nWF, 1);

gaussPts = [-1/sqrt(3), 1/sqrt(3)];
gaussWts = [1, 1];

for e = 1:fracture.nEdges

    xe = fracture.edge_centers(e,:)';
    Le = fracture.edge_lengths(e);

    nu = fracture.edge_normals(e,:)';
    nu = nu / norm(nu);

    % tangent direction along the edge
    t = cross([1;0;0], nu);
    t = t / norm(t);

    fluxIntegral = 0;

    for q = 1:2

        xq = xe + (Le/2) * gaussPts(q) * t;

        exact = exactFun(xq);

        fluxIntegral = fluxIntegral + gaussWts(q) * dot(exact.wFrac, nu);

    end

    wFExact(fracture.edge_flux_dofs(e)) = (Le/2) * fluxIntegral;
end

fprintf('\n[# MIXED-DIMENSIONAL FLUX ENERGY ERRORS]\n');

W  = flowBlocks.W;
WF = flowBlocks.WF;

dw  = w  - wExact;
dwF = wF - wFExact;

% (1) Bulk flux energy error
bulkFluxError = sqrt(dw' * W * dw);

% (2) Fracture flux energy error
fractureFluxError = sqrt(dwF' * WF * dwF);

fprintf('  Bulk flux energy error:      %.6e\n', bulkFluxError);

fprintf('  Fracture flux energy error:  %.6e\n', fractureFluxError);


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

    domegaDy = 2 * a * (1 - 2*y) * c^2;
    domegaDz = 2 * c * (1 - 2*z) * a^2;

    wx = -sign(x(1) - xFrac) * ((5/2) * d^(3/2) + omega);

    wy = -d * domegaDy;
    wz = -d * domegaDz;

    exact.wBulk = [wx; wy; wz];

    exact.wBulkMinus = [ omega; -d*domegaDy; -d*domegaDz ];
    exact.wBulkPlus  = [-omega; -d*domegaDy; -d*domegaDz ];

    exact.wFrac = [0; domegaDy; domegaDz];
end
