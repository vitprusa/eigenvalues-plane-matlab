clear all close all

% Isosceles right triangle

a = 0;
b = pi;
c = 0;
d = b; % Constraint for 45 degrees angle

x_range = [a, b];
y_range = [c, d];

% Grid spacing must be the same in both directions
% M is the number of DOF (interior grid points) in one slice (before applying the mask)

M = 120;
h = (b-a)/(M+1);

x_vec = x_range(1) : h : x_range(2);
y_vec = y_range(1) : h : y_range(2);

[X, Y] = ndgrid(x_vec, y_vec); 

% Defining function for the domain
phi = @(x,y) indicator_isosceles_triangle(x, y, a, b, c);   % Isosceles right triangle

% Create mask
global_mask = phi(X, Y) > 0; 

spy(global_mask)

%
% TESTING
%

% ANALYTICAL FORMULA
% We need a function that vanishes on boundary, and for which we can
% calculate the laplacian

mm = 1; % mm and nn must be different to have a non-trivial solution
nn = 3; 
u = @(x, y) sin(mm*(pi/b)*x).*sin(nn*(pi/b)*y) - sin(mm*(pi/b)*y).*sin(nn*(pi/b)*x) ;
laplace_u = @(x, y) -(mm^2 + nn^2)*(pi/b)^2*u(x, y);

u_grid = u(X, Y);
laplace_u_grid = laplace_u(X, Y);

u_masked_vals_grid = zeros(size(X));
u_masked_vals_grid(global_mask) = u_grid(global_mask);  

figure
surf(X,Y,u_masked_vals_grid)

laplace_u_masked_vals_grid = zeros(size(X));
laplace_u_masked_vals_grid(global_mask) = laplace_u_grid(global_mask); 

% figure
% surf(X,Y,laplace_u_masked_vals_grid)

% NUMERICAL COMPUTATION, DST
tic
laplace_u_masked_vals_grid_num = dst_laplace_grid(u_masked_vals_grid, h);

toc

% RESIDUAL
residual = laplace_u_masked_vals_grid - laplace_u_masked_vals_grid_num;
% Calculate the norm of the residual for convergence check
normResidual = norm(residual(global_mask));
disp(['Norm of the residual: ', num2str(normResidual)]);


