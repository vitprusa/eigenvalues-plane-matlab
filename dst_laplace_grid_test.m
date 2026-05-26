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

% Indexing, just for debugging
%seq = reshape(1:numel(X), size(X));

%masked_seq = NaN(size(X));
%masked_seq(global_mask) = seq(global_mask);
%masked_uvals_grid = masked_seq; 

%
% TESTING
%

% ANALYTICAL FORMULA
% We need a function that vanishes on boundary, and for which we can
% calculate the laplacian

mm = 5;
nn = 1;
u = @(x, y) sin(mm*x).*sin(nn*y);
laplace_u = @(x, y) (-mm.^2 - nn.^2)*u(x, y);

u_grid = u(X, Y);
laplace_u_grid = laplace_u(X, Y);

u_masked_vals_grid = NaN(size(X));
u_masked_vals_grid(global_mask) = u_grid(global_mask);  

laplace_u_masked_vals_grid = NaN(size(X));
laplace_u_masked_vals_grid(global_mask) = laplace_u_grid(global_mask);  

% NUMEIRCAL COMPUTATION, DST
tic
laplace_u_masked_vals_grid_num = dst_laplace_grid(u_masked_vals_grid, h);
toc

% RESIDUAL
residual = laplace_u_masked_vals_grid - laplace_u_masked_vals_grid_num;
% Calculate the norm of the residual for convergence check
normResidual = norm(residual(global_mask));
disp(['Norm of the residual: ', num2str(normResidual)]);
