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
M = 60;
h = (b-a)/(M+1);

x_vec = x_range(1) : h : x_range(2);
y_vec = y_range(1) : h : y_range(2);

% b = x_vec(3);
% c = x_vec(5);


[X, Y] = ndgrid(x_vec, y_vec); 

% Defining function for the domain
phi = @(x,y) indicator_ellipse(x, y, a, b, c, e);   % ellipse minus quadrant domain

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

% Eigenvalue solver, subset of eigenvalues, L in operator form
k = 10; % We want first k eigenvalues
%opts.issym  = true;           % operator is symmetric, BEWARE the matrix is not, but the operator is
%opts.isreal = true;           % set if A is real
%opts.tol = 1e-10;
%opts.maxit = 500;
%opts.SubspaceDimension = 100;
%opts.disp   = 1;

%'bothendsreal'
tic
D = eigs(Lop, dofs, k, 'largestreal', ...
    'Tolerance',           1e-10, ...
    'MaxIterations',       100, ...
    'SubspaceDimension',   800, ...
    'IsFunctionSymmetric', true, ...   
    'Display',             1);
%D = eigs(Lop, dofs, k, sigma, opts);
D = real(D);
toc
eigs_numerical_subset = sort(D, 'descend');

% % Assemble matrix whose i-th column is Lop applied to i-th column of
% % identity matrix in the degrees of freedom space. Can I vectorise this?
% idm = eye(dofs);
% L = zeros(dofs);
% parfor i = 1:dofs
%     L(:,i) = Lop(idm(:,i));
% end
% 
% % This does not work
% % L = Lop(eye(dofs));
% tic
% [~, D] = eig(L);
% eigs_numerical = sort(diag(D), 'descend');
% toc

% First 3 eigenvalues known in literature 

lambda = [-5.868746216295;
          -11.52599695049;
          -15.14021979035]; 

Lambda = [-5.87631;
	      -11.5438;
	      -15.1471;
	      -15.927;
	      -21.132;
	      -27.0968;
	      -27.5345;
	      -29.9545;
	      -34.7477;
	      -41.4342;
	      -42.2547;
	      -44.5974;
	      -47.0828;
	      -50.5084;
	      -52.4641;
	      -57.9654;
	      -61.6035;
	      -64.0446;
	      -66.2612;
	      -68.1086;
	      -73.9802;
	      -74.5986;
	      -78.2619;
	      -85.2798;
	      -85.8806;
	      -89.5776;
	      -91.5177;
	      -92.8236;
	      -99.254;
	      -101.106;
	      -102.616;
	      -105.597;
	      -112.505;
	      -113.355;
	      -115.266;
	      -116.428;
	      -121.837;
	      -123.354;
	      -130.005;
	      -130.168];
          
        
% Be careful, the linear ordering is not by n x n blocks!
% We rather match a block, not the first dof eigenvalues of the continuous
% operator!
% norm(eigs_numerical(1:40) - Lambda)
% norm(eigs_numerical(1:40) - Lambda, Inf)
norm(eigs_numerical_subset(1:k) - Lambda(1:k))
norm(eigs_numerical_subset(1:k) - Lambda(1:k), Inf)

