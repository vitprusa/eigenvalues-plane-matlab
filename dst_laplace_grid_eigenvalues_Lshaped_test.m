clc;
clear

% L-shaped domain 

a = -1;
b = 0;
c = 1;
d = -1;
e = 0;
f = 1;

x_range = [a, c];
y_range = [d, f];


% Grid spacing must be the same in both directions
% M is the number of DOF (interior grid points) in one slice (before applying the mask)
M = 23;
h = (c-a)/(M+1);

x_vec = x_range(1) : h : x_range(2);
y_vec = y_range(1) : h : y_range(2);
[X, Y] = ndgrid(x_vec, y_vec); 

% Defining function for the domain
phi = @(x,y) indicator_Lshape(x, y, a, b, c, d, e, f);   % L-shaped domain

% Create mask
global_mask = phi(X, Y) > 0; 

spy(global_mask)

% NUMERICAL COMPUTATION, DST
% Laplace operator acting on values vector
Lop = @(x) vals_grid_to_vals_vec( ...
    dst_laplace_grid( ...
    vals_vec_to_vals_grid(x, global_mask), ...
    global_mask, h), ...
    global_mask);

% Number of degrees of freedom
dofs = sum(global_mask, 'all');

% Assemble matrix whose i-th column is Lop applied to i-th column of
% identity matrix in the degrees of freedom space. Can I vectorise this?
idm = eye(dofs);
L = zeros(dofs);
for i = 1:dofs
    L(:,i) = Lop(idm(:,i));
end

% This does not work
% L = Lop(eye(dofs));

[~, D] = eig(L);
eigs_numerical = sort(diag(D), 'descend');

% FIRST TEN EIGENVALUES known in literature

lambda = [-9.6397238440;
          -15.1972519267;
          -19.7392088022; % 2*pi^2
          -29.5214811142;
          -31.9126359594;
          -41.4745098927;
          -44.9484877828;
          -49.3480220054; % 5*pi^2
          -49.3480220054; % 5*pi^2
          -56.7096098902];

% Be careful, the linear ordering is not by n x n blocks!
% We rather match a block, not the first dof eigenvalues of the continuous
% operator!
norm(eigs_numerical(1:10) - lambda)
norm(eigs_numerical(1:10) - lambda, Inf)
