clc;
clear

% Rectangular domain [0, pi] x [0, pi]
a = 0;
b = pi;
c = 0;
d = pi;

% Number of Chebyshev points in each direction
N = 20;
Nx = N;
Ny = N;

% Second derivative matrices
D2x = diffmat(Nx, 2, [a b]);
D2y = diffmat(Ny, 2, [c d]);

% Remove boundary points (first and last) for Dirichlet BCs
D2x_int = D2x(2:end-1, 2:end-1);
D2y_int = D2y(2:end-1, 2:end-1);

% 2D Laplacian via Kronecker product
Ix = speye(Nx - 2);
Iy = speye(Ny - 2);
L = kron(Iy, D2x_int) + kron(D2y_int, Ix);

fprintf('DOFs: %d\n', size(L, 1));

% Solve eigenvalue problem (eigenvalues of -Laplacian are positive)
[~, D] = eig(full(-L));
eigs_chebfun = sort(real(diag(D)), 'ascend');

% Analytical eigenvalues: (n*pi/(b-a))^2 + (m*pi/(d-c))^2 = n^2 + m^2
M = 50;
kk = 1:M;
eigs_all = kk.^2 + kk'.^2;
eigs_analytical = sort(eigs_all(:), 'ascend');

k = min(length(eigs_chebfun), length(eigs_analytical));
fprintf('Comparing first %d eigenvalues (Chebfun vs analytical):\n', k);
fprintf('  2-norm error: %e\n', norm(eigs_chebfun(1:k) - eigs_analytical(1:k)));
fprintf('  Inf-norm error: %e\n', norm(eigs_chebfun(1:k) - eigs_analytical(1:k), Inf));

figure
plot(1:length(eigs_chebfun), eigs_chebfun, 'o', 'DisplayName', 'Chebfun')
hold on
plot(1:k, eigs_analytical(1:k), 'x', 'DisplayName', 'Analytical')
hold off
xlabel('Eigenvalue index')
ylabel('Eigenvalue')
legend
title('Chebfun vs analytical eigenvalues (rectangle)')
