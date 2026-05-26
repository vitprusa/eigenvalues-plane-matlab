clc;
clear
addpath('..')

% Isosceles right triangle with vertices at (0,0), (pi,0), (pi,pi)
% Defined by: x > 0, x < pi, y > 0, y < x

% Define geometry using a polygon
a = 0;
b = pi;
pgon = [2; 3; a; b; b; a; a; b];  % triangle (a,a)-(b,a)-(b,b)

gd = pgon;
ns = char('T1')';
sf = 'T1';

[dl, bt] = decsg(gd, sf, ns);

% Create PDE model
model = createpde();
geometryFromEdges(model, dl);

% Boundary conditions: homogeneous Dirichlet
applyBoundaryCondition(model, 'dirichlet', 'Edge', 1:model.Geometry.NumEdges, 'u', 0);

% Coefficients for -Laplacian eigenvalue problem: -div(c*grad(u)) + a*u = lambda*d*u
specifyCoefficients(model, 'm', 0, 'd', 1, 'c', 1, 'a', 0, 'f', 0);

% Generate mesh
generateMesh(model, 'Hmax', 0.05, 'GeometricOrder', 'quadratic');

% Solve eigenvalue problem in range [0, 200]
result = solvepdeeig(model, [0, 200]);
eigs_fem = sort(result.Eigenvalues, 'ascend');

% Analytical eigenvalues for the isosceles right triangle with legs of length pi
% are a subset of -(n^2 + m^2) for n,m >= 1 (those with antisymmetric modes
% along the hypotenuse). Here we just compare with brute-force enumeration.
M = 50;
k = 1:M;
eigs_all = k.^2 + k'.^2;
eigs_analytical = sort(unique(eigs_all(:)), 'ascend');

k = min(length(eigs_fem), length(eigs_analytical));
fprintf('Comparing first %d eigenvalues (FEM vs analytical):\n', k);
fprintf('  2-norm error: %e\n', norm(eigs_fem(1:k) - eigs_analytical(1:k)));
fprintf('  Inf-norm error: %e\n', norm(eigs_fem(1:k) - eigs_analytical(1:k), Inf));
