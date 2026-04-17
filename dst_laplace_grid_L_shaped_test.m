clear all close all

% L-shaped domain 
a = 0;
c = pi;
d = 0;
f = pi;

x_range = [a, c];
y_range = [d, f];

% Grid spacing must be the same in both directions
% M is the number of DOF (interior grid points) in one slice (before applying the mask)
M = 43;
h = (c-a)/(M+1);

P = (M+1)/2;
Q = (M+1)/4;

b = P*h;
e = Q*h;

x_vec = x_range(1) : h : x_range(2);
y_vec = y_range(1) : h : y_range(2);

[X, Y] = ndgrid(x_vec, y_vec); 

% Defining function for the domain
phi = @(x,y) indicator_Lshape(x, y, a, b, c, d, e, f);   % L-shaped domain

% Create mask
global_mask = phi(X, Y) > 0; 

spy(global_mask)

%
% TESTING
%

% ANALYTICAL FORMULA
% We need a function that vanishes on boundary, and for which we can
% calculate the laplacian

mm = 2;
nn = 6;
u = @(x, y) sin(mm/(P*h)*pi*x).*sin(nn/(Q*h)*pi*y);
laplace_u = @(x, y) (-(pi).^2 / ((P*h).^2)*mm.^2 - (pi).^2 / ((Q*h).^2)*nn.^2)*u(x, y);


u_grid = u(X, Y);
laplace_u_grid = laplace_u(X, Y);

u_masked_vals_grid = zeros(size(X));
u_masked_vals_grid(global_mask) = u_grid(global_mask);  

figure
surf(X,Y,u_masked_vals_grid)

laplace_u_masked_vals_grid = zeros(size(X));
laplace_u_masked_vals_grid(global_mask) = laplace_u_grid(global_mask); 

figure
surf(X,Y,laplace_u_masked_vals_grid)

% NUMERICAL COMPUTATION, DST
tic
laplace_u_masked_vals_grid_num = dst_laplace_grid(u_masked_vals_grid, global_mask, h);

toc

% RESIDUAL
residual = laplace_u_masked_vals_grid - laplace_u_masked_vals_grid_num;
% Calculate the norm of the residual for convergence check
normResidual = norm(residual(global_mask));
disp(['Norm of the residual: ', num2str(normResidual)]);
