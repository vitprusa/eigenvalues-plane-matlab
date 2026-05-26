clc;
clear

% H domain bounding box: [-1, 2] x [-2, 1]
x_range = [-1, 2];
y_range = [-2, 1];

% Grid spacing must be the same in both directions
% M is the number of DOF (interior grid points) in one slice (before applying the mask)
M = 71;
h = (x_range(2) - x_range(1))/(M+1);

x_vec = x_range(1) : h : x_range(2);
y_vec = y_range(1) : h : y_range(2);
[X, Y] = ndgrid(x_vec, y_vec);

% Defining function for the domain
phi = @(x,y) indicator_H(x, y);

% Create mask
global_mask = phi(X, Y) > 0;

%spy(global_mask)

% NUMERICAL COMPUTATION, DST
% Laplace operator acting on values vector
Lop = @(x) vals_grid_to_vals_vec( ...
    dst_laplace_grid( ...
    vals_vec_to_vals_grid(x, global_mask), h), ...
    global_mask);

% Number of degrees of freedom
dofs = sum(global_mask, 'all');

% Assemble matrix whose i-th column is Lop applied to i-th column of
% identity matrix in the degrees of freedom space.
idm = eye(dofs);
L = zeros(dofs);
parfor i = 1:dofs
    L(:,i) = Lop(idm(:,i));
end

[~, D] = eig(L);
eigs_numerical = sort(diag(D), 'descend');
