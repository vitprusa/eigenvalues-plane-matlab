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

% Solve eigenvalue problem in range [0, 200]
result = solvepdeeig(model, [0, 200]);
eigs_fem = sort(result.Eigenvalues, 'ascend');

% Analytical eigenvalues: (n*pi/(b-a))^2 + (m*pi/(d-c))^2 = n^2 + m^2
% for the [0,pi] x [0,pi] square
M = 50;
k = 1:M;
eigs_all = k.^2 + k'.^2;
eigs_analytical = sort(eigs_all(:), 'ascend');

k = min(length(eigs_fem), length(eigs_analytical));
fprintf('Comparing first %d eigenvalues (FEM vs analytical):\n', k);
fprintf('  2-norm error: %e\n', norm(eigs_fem(1:k) - eigs_analytical(1:k)));
fprintf('  Inf-norm error: %e\n', norm(eigs_fem(1:k) - eigs_analytical(1:k), Inf));
