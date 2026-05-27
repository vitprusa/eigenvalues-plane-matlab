clear all close all

% Resolution of the BVP on L-shaped domain

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

P = (M+1)/2; % Must be factor of M+1
Q = (M+1)/4; % Must be factor of M+1

b = P*h;
e = Q*h;

x_vec = x_range(1) : h : x_range(2);
y_vec = y_range(1) : h : y_range(2);

[X, Y] = ndgrid(x_vec, y_vec); 

% Defining function for the domain
phi = @(x,y) indicator_Lshape(x, y, a, b, c, d, e, f);   % L-shaped domain

% Create mask
global_mask = phi(X, Y) > 0; 

%spy(global_mask)



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


% ANALYTICAL FORMULA
% We need a test function that vanishes on boundary, and for which we can
% calculate the laplacian

% TEST FUNCTIONS

% 1) Eigenfunction

% mm = 2; % Must be less or equal to P
% nn = 3; % Must be less or equal to Q
% u = @(x, y) sin(mm/(P*h)*pi*x).*sin(nn/(Q*h)*pi*y);
% laplace_u = @(x, y) (-(pi).^2 / ((P*h).^2)*mm.^2 - (pi).^2 / ((Q*h).^2)*nn.^2)*u(x, y);

% 2) Polynomial

u = @(x, y) x.*(x-pi).*(x-P*h).*y.*(y-pi).*(y-Q*h);
laplace_u = @(x, y) 2*(3*x-pi-P*h).*y.*(y-pi).*(y-Q*h)+2*(3*y-pi-Q*h).*x.*(x-pi).*(x-P*h);

u_grid = u(X, Y);
laplace_u_grid = laplace_u(X, Y);

u_masked_vals_grid = zeros(size(X));
u_masked_vals_grid(global_mask) = u_grid(global_mask); 


% Solve the boundary value problem: L u_vec = b_vec
b_vec = laplace_u_grid(global_mask);
u_vec = L\b_vec;
% Reshape the solution vector back to the grid
u_sol_grid = zeros(size(X));
u_sol_grid(global_mask) = u_vec;

% Calculate the error between the numerical solution and the analytical solution
error = u_sol_grid - u_masked_vals_grid;

% Calculate the maximum error and display it
maxError = max(abs(error(global_mask)));
disp(['Maximum error: ', num2str(maxError)]);

figure
surf(X,Y,u_sol_grid)
colorbar;
xlabel('$x$','Interpreter','latex');
ylabel('$y$','Interpreter','latex');
zlabel('$u(x,y)$','Interpreter','latex');
title('Numerical Solution on L-shaped Domain');

