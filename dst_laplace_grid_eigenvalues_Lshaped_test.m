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
M = 60;
h = (b-a)/(M+1);

x_vec = x_range(1) : h : x_range(2);
y_vec = y_range(1) : h : y_range(2);
[X, Y] = ndgrid(x_vec, y_vec); 

% Defining function for the domain
phi = @(x,y) indicator_Lshape(x, y, a, b, c, d, e, f);   % L-shaped domain

% Create mask
tic
global_mask = phi(X, Y) > 0; 
toc

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
% identity matrix in the degrees of freedom space. Can I vectorise this?
idm = eye(dofs);
L = zeros(dofs);
parfor i = 1:dofs
    L(:,i) = Lop(idm(:,i));
end

% This does not work
% L = Lop(eye(dofs));

% CHECK: With some large M, L is no longer perfectly symmetric beacuse of 
% numerical errors

% norm(L-L','fro')


%DD = eigs(Lop,dofs,10,'smallestabs');
% [~, D] = eig(L);
tic
D = eig(L);
toc
%eigs_numerical = sort(diag(D), 'descend');
%eigs_numerical = sort(D,'descend');

% We have to consider the real part of D, otherwise the sorting does not
% work well
eigs_numerical = sort(real(D),'descend'); 


% First 10 eigenvalues known in literature

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

% First 40 eigenvalues

Lambda = [-9.65934;
	      -15.1988;
	      -19.7414;
	      -29.5285;
	      -31.967;
	      -41.526;
	      -44.9729;
	      -49.3775;
	      -49.38;
	      -56.7945;
	      -65.4565;
	      -71.1954;
	      -71.6685;
	      -79.0801;
	      -89.5956;
	      -92.4954;
	      -97.6157;
	      -98.9282;
	      -98.9364;
	      -101.916;
	      -112.734;
	      -115.911;
	      -128.783;
	      -128.834;
	      -130.71;
	      -130.894;
	      -143.209;
	      -151.945;
	      -155.561;
	      -163.203;
	      -165.685;
	      -165.901;
	      -168.777;
	      -168.895;
	      -178.933;
	      -181.038;
	      -185.429;
	      -199.12;
	      -199.158;
	      -203.332];

% Be careful, the linear ordering is not by n x n blocks!
% We rather match a block, not the first dof eigenvalues of the continuous
% operator!
norm(eigs_numerical(1:40) - Lambda)
norm(eigs_numerical(1:40) - Lambda, Inf)
