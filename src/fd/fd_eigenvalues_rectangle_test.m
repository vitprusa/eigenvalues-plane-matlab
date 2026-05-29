clc;
clear

% Rectangular domain [0, pi] x [0, pi]
a = 0;
b = pi;
c = 0;
d = pi;

% Number of interior grid points in each direction
M = 40;
hx = (b - a) / (M + 1);
hy = (d - c) / (M + 1);

% 1D second derivative matrix (interior points only, Dirichlet BCs)
e = ones(M, 1);
D2x = spdiags([e, -2*e, e], [-1, 0, 1], M, M) / hx^2;
D2y = spdiags([e, -2*e, e], [-1, 0, 1], M, M) / hy^2;

% 2D Laplacian via Kronecker product: I ⊗ D2x + D2y ⊗ I
Ix = speye(M);
Iy = speye(M);
L = kron(Iy, D2x) + kron(D2y, Ix);

fprintf('DOFs: %d\n', size(L, 1));

% Solve eigenvalue problem (eigenvalues of -Laplacian are positive)
[~, D] = eig(full(-L));
eigs_fd = sort(diag(D), 'ascend');

% Analytical eigenvalues: (n*pi/(b-a))^2 + (m*pi/(d-c))^2 = n^2 + m^2
k = 1:M^2;
eigs_all = k.^2 + k'.^2;
eigs_analytical = sort(eigs_all(:), 'ascend');

k = min(length(eigs_fd), length(eigs_analytical));
fprintf('Comparing first %d eigenvalues (FD vs analytical):\n', k);
fprintf('  2-norm error: %e\n', norm(eigs_fd(1:k) - eigs_analytical(1:k)));
fprintf('  Inf-norm error: %e\n', norm(eigs_fd(1:k) - eigs_analytical(1:k), Inf));

figure
plot(1:length(eigs_fd), eigs_fd, 'o', 'DisplayName', 'FD')
hold on
plot(1:k, eigs_analytical(1:k), 'x', 'DisplayName', 'Analytical')
hold off
xlabel('Eigenvalue index')
ylabel('Eigenvalue')
legend
title('FD vs analytical eigenvalues (rectangle)')
