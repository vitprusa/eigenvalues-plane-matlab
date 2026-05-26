clc;
clear

% H-shaped domain composed of three rectangles:
%   left column:  [-1, 0] x [-2, 1]
%   center bar:   [ 0, 1] x [-1, 0]
%   right column: [ 1, 2] x [-2, 1]

% Define the three rectangles as columns [x_left; x_right; y_bottom; y_top]
R1 = [3; 4; -1; 0; 0; 1; -2; -2; 1; 1];   % left column
R2 = [3; 4;  0; 1; 1; 0; -1; -1; 0; 0];   % center bar
R3 = [3; 4;  1; 2; 2; 1; -2; -2; 1; 1];   % right column

% Pad to equal length and combine
gd = [R1, R2, R3];
ns = char('R1', 'R2', 'R3')';
sf = 'R1 + R2 + R3';

[dl, bt] = decsg(gd, sf, ns);

% Create PDE model
model = createpde();
geometryFromEdges(model, dl);

% Boundary conditions: homogeneous Dirichlet
applyBoundaryCondition(model, 'dirichlet', 'Edge', 1:model.Geometry.NumEdges, 'u', 0);

% Coefficients for -Laplacian eigenvalue problem: -div(c*grad(u)) + a*u = lambda*d*u
specifyCoefficients(model, 'm', 0, 'd', 1, 'c', 1, 'a', 0, 'f', 0);

% Generate mesh
generateMesh(model, 'Hmax', 0.05, 'GeometricOrder', 'quadratic');

% Solve eigenvalue problem in range [0, 200]
result = solvepdeeig(model, [0, 200]);
eigs_fem = sort(result.Eigenvalues, 'ascend');

% First 40 known eigenvalues (Wolfram Language, positive convention)
Lambda_known = [7.77338;
    8.58911;
    13.9507;
    13.9546;
    14.3126;
    17.7206;
    19.748;
    24.8349;
    26.4306;
    26.4816;
    33.0826;
    37.2337;
    37.303;
    39.2822;
    40.8027;
    42.3474;
    45.3896;
    46.4166;
    47.6086;
    49.4738;
    49.4815;
    52.1734;
    55.1058;
    58.9657;
    62.9425;
    63.0716;
    63.4004;
    65.7973;
    67.8505;
    72.4734;
    77.4846;
    79.4798;
    80.3673;
    80.4442;
    85.8397;
    86.6899;
    88.6362;
    92.0113;
    92.9729;
    96.3744];

k = min(length(eigs_fem), length(Lambda_known));
fprintf('Comparing first %d eigenvalues:\n', k);
fprintf('  2-norm error: %e\n', norm(eigs_fem(1:k) - Lambda_known(1:k)));
fprintf('  Inf-norm error: %e\n', norm(eigs_fem(1:k) - Lambda_known(1:k), Inf));
