clc;
clear

% H domain bounding box: [-1, 2] x [-2, 1]
x_range = [-1, 2];
y_range = [-2, 1];

% Grid spacing must be the same in both directions
% M is the number of DOF (interior grid points) in one slice (before applying the mask)
M = 68;
h = (x_range(2) - x_range(1))/(M+1);

x_vec = x_range(1) : h : x_range(2);
y_vec = y_range(1) : h : y_range(2);
[X, Y] = ndgrid(x_vec, y_vec);

% Defining function for the domain
phi = @(x,y) indicator_H(x, y);

% Create mask
global_mask = phi(X, Y) > 0;

spy(global_mask)

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

% First 40 eigenvalues, Wolfram Language
Lambda = -[7.77338;
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

% Be careful, the linear ordering is not by n x n blocks!
% We rather match a block, not the first dof eigenvalues of the continuous
% operator!
norm(eigs_numerical(1:40) - Lambda)
norm(eigs_numerical(1:40) - Lambda, Inf)