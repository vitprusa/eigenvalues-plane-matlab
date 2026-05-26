clc;
clear

% Domain [a, b] x [c, d]
a = 0;
b = pi;
c = 0;
d = pi;

x_range = [a, b];
y_range = [c, d];

% Grid spacing must be the same in both directions
% M is the number of DOF (interior grid points) in one slice (before applying the mask)
M = 20;
h = (b-a)/(M+1);

x_vec = x_range(1) : h : x_range(2);
y_vec = y_range(1) : h : y_range(2);
[X, Y] = ndgrid(x_vec, y_vec); 

% Defining function for the domain
phi = @(x,y) indicator_rectangle(x, y, a, b, c, d);   % rectangle

% Create mask
global_mask = phi(X, Y) > 0; 


% NUMERICAL COMPUTATION, DST
% Laplace operator acting on values vector
Lop = @(x) vals_grid_to_vals_vec( ...
    dst_laplace_grid( ...
    vals_vec_to_vals_grid(x, global_mask), h), ...
    global_mask);

% Number of degrees of freedom
dofs = sum(global_mask, 'all');

% Assemble matrix whose i-th column is Lop applied to i-th column of
% identity matrix in the degeers of freedom space. Can I vectorise this?
idm = eye(dofs);
L = zeros(dofs);
for i = 1:dofs
    L(:,i) = Lop(idm(:,i));
end

% This does not work
% L = Lop(eye(dofs));

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
