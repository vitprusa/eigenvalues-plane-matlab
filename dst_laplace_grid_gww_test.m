clear all close all

% GWW1 and GWW2 domain 

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
M = 11;
h = (b-a)/(M+1);

x_vec = x_range(1) : h : x_range(2);
y_vec = y_range(1) : h : y_range(2);

% b = x_vec(3);
% c = x_vec(5);


[X, Y] = ndgrid(x_vec, y_vec); 

% Defining function for the domain: choose which one between GWW1 and GWW2
%phi = @(x,y) indicator_gww1(x, y, a, b, c, d, f, g);     % GWW1 domain
phi = @(x,y) indicator_gww2(x, y, a, b, c, d, f, g, i);   % GWW2 domain

% Create mask
global_mask = phi(X, Y) > 0; 

spy(global_mask)

%
% TESTING
%

% ANALYTICAL FORMULA
% We need a function that vanishes on boundary, and for which we can
% calculate the laplacian

% mm = 1; 
% nn = 1; 
% u = @(x, y) sin(mm*(pi/h)*x).*sin(nn*(pi/h)*y);
% laplace_u = @(x, y) (-mm.^2 - nn.^2)*(pi/h)^2*u(x, y);
% 
% 
% u_grid = u(X, Y);
% laplace_u_grid = laplace_u(X, Y);
% 
% u_masked_vals_grid = zeros(size(X));
% u_masked_vals_grid(global_mask) = u_grid(global_mask);  
% 
% figure
% surf(X,Y,u_masked_vals_grid)
% 
% laplace_u_masked_vals_grid = zeros(size(X));
% laplace_u_masked_vals_grid(global_mask) = laplace_u_grid(global_mask); 
% 
% figure
% surf(X,Y,laplace_u_masked_vals_grid)
% 
% % NUMERICAL COMPUTATION, DST
% tic
% laplace_u_masked_vals_grid_num = dst_laplace_grid(u_masked_vals_grid, h);
% 
% toc
% 
% % RESIDUAL
% residual = laplace_u_masked_vals_grid - laplace_u_masked_vals_grid_num;
% % Calculate the norm of the residual for convergence check
% normResidual = norm(residual(global_mask));
% disp(['Norm of the residual: ', num2str(normResidual)]);
