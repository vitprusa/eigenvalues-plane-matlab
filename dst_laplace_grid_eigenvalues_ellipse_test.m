clear all close all

% Three quarters of an ellipse domain 

a = -2;
b = 0;
c = 2;
d = -1;
e = 0;
f = 1;

x_range = [a, c];
y_range = [d, f];

% Grid spacing must be the same in both directions
% M is the number of DOF (interior grid points) in one slice (before applying the mask)
M = 25;
h = (b-a)/(M+1);

x_vec = x_range(1) : h : x_range(2);
y_vec = y_range(1) : h : y_range(2);

% b = x_vec(3);
% c = x_vec(5);


[X, Y] = ndgrid(x_vec, y_vec); 

% Defining function for the domain
phi = @(x,y) indicator_ellipse(x, y, a, b, c, e);   % 3/4 ellipse domain

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
parfor i = 1:dofs
    L(:,i) = Lop(idm(:,i));
end

% This does not work
% L = Lop(eye(dofs));
tic
[~, D] = eig(L);
eigs_numerical = sort(diag(D), 'descend');
toc

% First 3 eigenvalues known in literature 

lambda = [-5.868746216295;
          -11.52599695049;
          -15.14021979035]; 
          
        
% Be careful, the linear ordering is not by n x n blocks!
% We rather match a block, not the first dof eigenvalues of the continuous
% operator!
norm(eigs_numerical(1:3) - lambda)
norm(eigs_numerical(1:3) - lambda, Inf)


