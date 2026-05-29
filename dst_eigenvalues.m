clc;
clear

% Domain [a, b] x [c, d]
a = 0;
b = pi;
c = 0;
d = pi;

x_range = [a, b];
y_range = [c, d];

% M is the number of DOF (interior grid points) in one slice (before applying the mask)
M = 40;

% Defining function for the domain
phi = @(x,y) indicator_rectangle(x, y, a, b, c, d);   % rectangle

% NUMERICAL COMPUTATION, DST
% Assemble the masked DST-based Laplace matrix and the grid/domain info
% [L, info] = make_dst_laplace_mat(x_range, y_range, M, phi);
[L, info] = make_dst_laplace_mat_batched(x_range, y_range, M, phi);

[~, D] = eig(L);
eigs_numerical = sort(diag(D), 'descend');

% ANALYTICAL FORMULA
% eigenvalues are given as - (n^2 +m^2)
k = 1:M; % This is overkill, but who cares
eigs_analytical = -(k.^2 + k'.^2);
eigs_analytical = sort(eigs_analytical(:), 'descend');

% Be careful, the linear ordering is not by n x n blocks!
% We rather match a block, not the first dof eigenvalues of the continucous
% operator!
norm(eigs_numerical - eigs_analytical)
norm(eigs_numerical - eigs_analytical, Inf)
