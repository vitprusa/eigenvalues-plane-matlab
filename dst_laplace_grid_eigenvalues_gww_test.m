clear all close all

% GWW1 and GWW2 domain (same eigenvalues)

a = -3;
b = -1;
c = 1;
d = 3;
e = -3;
f = -1;
g = 1;
i = 3;

x_range = [a, d];
y_range = [e, i];

% Grid spacing must be the same in both directions
% M is the number of DOF (interior grid points) in one slice (before applying the mask)
M = 25;
h = (b-a)/(M+1);

x_vec = x_range(1) : h : x_range(2);
y_vec = y_range(1) : h : y_range(2);

% b = x_vec(3);
% c = x_vec(5);


[X, Y] = ndgrid(x_vec, y_vec); 

% Defining function for the domain: choose which one between GWW1 and GWW2
%phi = @(x,y) indicator_global_gww1(x, y, a, b, c, d, f, g);     % GWW1 domain
phi = @(x,y) indicator_global_gww2(x, y, a, b, c, d, f, g, i);   % GWW2 domain

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

% First 25 eigenvalues known in literature 
% (same eigenvalues for GWW1 and GWW2)

lambda = [-2.53794399980;
          -3.65550971352;
          -5.17555935622; 
          -6.53755744376;
          -7.24807786256;
          -9.20929499840;
          -10.5969856913;
          -11.5413953956; 
          -12.3370055014; 
          -13.0536540557;
          -14.3138624643;
          -15.8713026200;
          -16.9417516880;
          -17.6651184368;
          -18.9810673877;
          -20.8823950433;
          -21.2480051774;
          -22.2328517930;
          -23.7112974848;
          -24.4792340693;
          -24.6740110027;
          -26.0802400997;
          -27.3040189211;
          -28.1751285815;
          -29.5697729132];


% Be careful, the linear ordering is not by n x n blocks!
% We rather match a block, not the first dof eigenvalues of the continuous
% operator!
norm(eigs_numerical(1:25) - lambda)
norm(eigs_numerical(1:25) - lambda, Inf)


