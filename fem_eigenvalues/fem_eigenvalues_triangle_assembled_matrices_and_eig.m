clc;
clear

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

% Analytical eigenvalues for the isosceles right triangle with legs of length pi
% are a subset of -(n^2 + m^2) for n,m >= 1 (those with antisymmetric modes
% along the hypotenuse). Here we just compare with brute-force enumeration.
% TODO
M = 90;
k = 1:M;
eigs_all = k.^2 + k'.^2;
eigs_analytical = sort(unique(eigs_all(:)), 'ascend');

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
title('FEM vs analytical eigenvalues (triangle)')
