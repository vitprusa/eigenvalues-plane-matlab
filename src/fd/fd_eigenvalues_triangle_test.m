clc;
clear

% Isosceles right triangle with vertices at (0,0), (pi,0), (pi,pi)
% Defined by: x > 0, x < pi, y > 0, y < x
a = 0;
b = pi;

x_range = [a, b];
y_range = [a, b];

% Number of interior grid points per side
M = 40;
h = (b - a) / (M + 1);

x_vec = x_range(1) : h : x_range(2);
y_vec = y_range(1) : h : y_range(2);
[X, Y] = ndgrid(x_vec, y_vec);

% Mask for interior points of the triangle: y < x
phi = @(x, y) indicator_isosceles_triangle(x, y, a, b, a);
global_mask = phi(X, Y) > 0;

dofs = sum(global_mask, 'all');
fprintf('DOFs: %d\n', dofs);

% Map between grid indices and DOF indices
dof_map = zeros(size(global_mask));
dof_map(global_mask) = 1:dofs;

% Assemble 2D Laplacian with 5-point stencil
L = sparse(dofs, dofs);
[m, n] = size(global_mask);
for j = 1:n
    for i = 1:m
        if ~global_mask(i, j)
            continue
        end
        idx = dof_map(i, j);
        L(idx, idx) = -4 / h^2;
        if i > 1 && global_mask(i-1, j)
            L(idx, dof_map(i-1, j)) = 1 / h^2;
        end
        if i < m && global_mask(i+1, j)
            L(idx, dof_map(i+1, j)) = 1 / h^2;
        end
        if j > 1 && global_mask(i, j-1)
            L(idx, dof_map(i, j-1)) = 1 / h^2;
        end
        if j < n && global_mask(i, j+1)
            L(idx, dof_map(i, j+1)) = 1 / h^2;
        end
    end
end

% Solve eigenvalue problem (eigenvalues of -Laplacian are positive)
[~, D] = eig(full(-L));
eigs_fd = sort(diag(D), 'ascend');

% Analytical eigenvalues for the isosceles right triangle
% (subset of n^2 + m^2 for the square)
k = 1:M^2;
eigs_all = k.^2 + k'.^2;
eigs_analytical = sort(unique(eigs_all(:)), 'ascend');

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
title('FD vs analytical eigenvalues (triangle)')
