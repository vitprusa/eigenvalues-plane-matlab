clc;
clear

% Domain [a, b] x [c, d]
a = 0;
b = pi;
c = 0; 
d = pi;

% Grid spacing must be the same in both directions
% M is the number of DOF in the x-direction
% Number of gridpoints (including end points) is M+2, grid spacing is h
M = 10000;
h = (b-a)/(M+1);

% Now we create the grid, including the endpoints
xfull = a:h:b; 
yfull = c:h:d;

x = xfull(2:end-1);
y = yfull(2:end-1);

% This is our grid, inner points only
[X,Y] = ndgrid(x,y);

% Differentiation on slice
kvec = @(x) -((1:length(x)).^2);
D2op = @(x) idst(kvec(x)' .* dst(x));

m = 5;
n = 6;
f = @(x, y) sin(m*x).*sin(n*y);
F = f(X, Y);

firstColumn = F(:, 1);
D2firstColumn = D2op(firstColumn);

firstRow = F(1, :);
D2firstRow = D2op(firstRow');

Lop = @(x) dst_laplacian_rectangle(x, X, Y, a, b, c, d);

laplaceToFvec = Lop(F(:));
reshape(laplaceToFvec, size(X));

% Compute the residual of the Laplacian operator applied to F
residual = reshape(laplaceToFvec + (m.^2 + n.^2) * F(:), size(X));
