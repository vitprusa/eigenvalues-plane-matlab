clc;
clear

% H domain bounding box: [-1, 2] x [-2, 1]
x_range = [-1, 2];
y_range = [-2, 1];

% Number of interior grid points per side (in x direction)
M = 40;
h = (x_range(2) - x_range(1)) / (M + 1);

x_vec = x_range(1) : h : x_range(2);
y_vec = y_range(1) : h : y_range(2);
[X, Y] = ndgrid(x_vec, y_vec);

% Mask for interior points of the H domain
phi = @(x, y) indicator_H(x, y);
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

% First 40 known eigenvalues (Wolfram Language, positive convention)
Lambda_known = [7.77338;
    8.58911;
    13.9507;
    13.9546;
    14.3126;
    17.7206;
    19.748;
    24.8349;
    26.4306;
    26.4816;
    33.0826;
    37.2337;
    37.303;
    39.2822;
    40.8027;
    42.3474;
    45.3896;
    46.4166;
    47.6086;
    49.4738;
    49.4815;
    52.1734;
    55.1058;
    58.9657;
    62.9425;
    63.0716;
    63.4004;
    65.7973;
    67.8505;
    72.4734;
    77.4846;
    79.4798;
    80.3673;
    80.4442;
    85.8397;
    86.6899;
    88.6362;
    92.0113;
    92.9729;
    96.3744];

k = min(length(eigs_fd), length(Lambda_known));
fprintf('Comparing first %d eigenvalues:\n', k);
fprintf('  2-norm error: %e\n', norm(eigs_fd(1:k) - Lambda_known(1:k)));
fprintf('  Inf-norm error: %e\n', norm(eigs_fd(1:k) - Lambda_known(1:k), Inf));

figure
plot(1:length(eigs_fd), eigs_fd, 'o', 'DisplayName', 'FD')
hold on
plot(1:k, Lambda_known(1:k), 'x', 'DisplayName', 'Wolfram')
hold off
xlabel('Eigenvalue index')
ylabel('Eigenvalue')
legend
title('FD vs known eigenvalues (H domain)')
