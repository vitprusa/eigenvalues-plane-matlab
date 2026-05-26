clc;
clear
addpath('..')

% Rectangular domain [0, pi] x [0, pi]
a = 0;
b = pi;
c = 0;
d = pi;

% Define geometry as a rectangle [x1 x2 x3 x4; y1 y2 y3 y4]
pgon = [3; 4; a; b; b; a; c; c; d; d];

gd = pgon;
ns = char('R1')';
sf = 'R1';

[dl, bt] = decsg(gd, sf, ns);

% Create PDE model
model = createpde();
geometryFromEdges(model, dl);

% Boundary conditions: homogeneous Dirichlet
applyBoundaryCondition(model, 'dirichlet', 'Edge', 1:model.Geometry.NumEdges, 'u', 0);

% Coefficients for -Laplacian eigenvalue problem: -div(c*grad(u)) + a*u = lambda*d*u
specifyCoefficients(model, 'm', 0, 'd', 1, 'c', 1, 'a', 0, 'f', 0);

% Generate mesh
generateMesh(model, 'Hmax', 0.1, 'GeometricOrder', 'quadratic');

% Assemble stiffness and mass matrices
FEM_raw = assembleFEMatrices(model, 'KM');
FEM_ns  = assembleFEMatrices(model, 'nullspace');
B = FEM_ns.B;
K = B' * FEM_raw.K * B;
M = B' * FEM_raw.M * B;

fprintf('DOFs: %d\n', size(K, 1));

% Solve generalized eigenvalue problem K*u = lambda*M*u
[~, D] = eig(full(K), full(M));
eigs_fem = sort(diag(D), 'ascend');

% Analytical eigenvalues: (n*pi/(b-a))^2 + (m*pi/(d-c))^2 = n^2 + m^2
% for the [0,pi] x [0,pi] square
M = 70;
k = 1:M;
eigs_all = k.^2 + k'.^2;
eigs_analytical = sort(eigs_all(:), 'ascend');

k = min(length(eigs_fem), length(eigs_analytical));
fprintf('Comparing first %d eigenvalues (FEM vs analytical):\n', k);
fprintf('  2-norm error: %e\n', norm(eigs_fem(1:k) - eigs_analytical(1:k)));
fprintf('  Inf-norm error: %e\n', norm(eigs_fem(1:k) - eigs_analytical(1:k), Inf));

figure
plot(1:length(eigs_fem), eigs_fem, 'o', 'DisplayName', 'FEM')
hold on
plot(1:k, eigs_analytical(1:k), 'x', 'DisplayName', 'Analytical')
hold off
xlabel('Eigenvalue index')
ylabel('Eigenvalue')
legend
title('FEM vs analytical eigenvalues (rectangle)')
