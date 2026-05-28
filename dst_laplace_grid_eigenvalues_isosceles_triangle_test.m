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

M = 40;
h = (b-a)/(M+1);


x_vec = x_range(1) : h : x_range(2);
y_vec = y_range(1) : h : y_range(2);

[X, Y] = ndgrid(x_vec, y_vec); 


% Defining function for the domain
phi = @(x,y) indicator_isosceles_triangle(x, y, a, b, c);   % Isosceles right triangle

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
for i = 1:dofs
    L(:,i) = Lop(idm(:,i));
end

% This does not work
% L = Lop(eye(dofs));
tic
% [~, D] = eig(L);
% eigs_numerical = sort(diag(D), 'descend');
D = eig(L);
eigs_numerical = sort(real(D), 'descend');
toc

% ANALYTICAL FORMULA
% eigenvalues are given as - (pi/b)^2*(n^2 + m^2)
k = 1:M; % This is overkill, but who cares
eigs_analytical_full = zeros(M+2);
eigs_analytical_full(2:(end-1),2:(end-1)) = -(pi/b)^2*(k.^2+k'.^2);
eigs_analytical_new = zeros(M+2);
eigs_analytical_new(global_mask) = eigs_analytical_full(global_mask);
eigs_analytical_new = sort(eigs_analytical_new(:), 'descend');
eigs_analytical_final = eigs_analytical_new(end-dofs+1:end);


% Be careful, the linear ordering is not by n x n blocks!
% We rather match a block, not the first dof eigenvalues of the continuous
% operator!
% norm(eigs_numerical - eigs_analytical_final)
norm(eigs_numerical - eigs_analytical_final, Inf)

% If we want to display the plot with both eigs_analytical_full and
% eigs_numerical 


% data_analytical_full = -(k.^2+k'.^2);
% data_numerical = eigs_numerical';
% data_analytical_full = sort(data_analytical_full(:),'descend');
% data_x1 = 1:numel(data_analytical_full);
% data_x2 = 1:numel(data_numerical);
% plot(data_analytical_full,ones(size(data_analytical_full)),'o','LineStyle','none')
% hold on
% plot(data_numerical,ones(size(data_numerical)),'o','LineStyle','none','MarkerEdgeColor','r')
% hold off


% Asymptotic behavior:

Area = (b-a)*(d-c)/2;
%kmax = round(0.8*dofs);
n = (1 : dofs)';
% n = (dofs - k +1: dofs)';
E = -eigs_numerical./n;
figure
plot(n,E,'o','LineStyle','none','MarkerEdgeColor','r')
%plot(n,E,'r')
hold on
plot(n,ones(1,dofs)*4*pi/Area,'k')
hold on

C = 4*pi/Area;              % asymptotic value
diffE = abs(E - C);         % absolute difference
% Tolerance: relative to median(E) or an absolute small number
tol = 0.047 * median(abs(E));   % it's very sensitive and it has to be adjusted
% Find first index where diff exceeds tolerance
% The control starts at sufficiently large value of n, 
% as we are interested in the high-order indexes 
n_c = find(diffE(0.6*dofs:end) > tol, 1, 'first');
n_critical = 0.6*dofs + n_c;

if isempty(n_critical)
    fprintf('No deviation found above tol = %g\n', tol)
else
    fprintf('Deviation starts at n = %d (tol = %g)\n', n_critical, tol)
    xline(n_critical, '--k'); 
end