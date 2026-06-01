% DST_LAPLACE_PARTIAL_SPECTRUM_SMALL_RECTANGLE First few Dirichlet-Laplacian eigenvalues on a small rectangle.
%
%   Computes only the first few eigenvalues of the discrete Laplace operator
%   with homogeneous Dirichlet boundary conditions on a small rectangle
%   [0, pi/2] x [0, pi/4] sitting inside the bounding box [0, pi] x [0, pi].
%   The domain is described by an indicator function embedded in a rectangular
%   bounding box and discretised on a uniform grid; the DST-based Laplace
%   operator is built with make_dst_laplace_op_batched and the leading
%   eigenvalues are obtained with the iterative solver eigs.
%
%   Generated from the dst_laplace_partial_spectrum template.

clc;
clear

% This script lives in results/; put the project sources on the path.
run(fullfile(fileparts(fileparts(mfilename('fullpath'))), 'startup.m'));

% Bounding box dimensions [a, b] x [c, d]
a = 0;
b = pi;
c = 0;
d = pi;

[x_range, y_range] = bounding_box(a, b, c, d);

% Number of DOF (interior grid points) in one slice, full bounding box.
% M+1 is a multiple of 4 so the rectangle edges pi/2 and pi/4 fall on grid lines.
M = 199;

% Number of eigenvalues (first few) to compute
k = 10;

% Dimension of the Krylov subspace used by eigs (tuning option). Must satisfy
% k < subspace_dim <= dofs; larger values can improve convergence on harder
% domains at the cost of more memory and work per iteration.
subspace_dim = 100;

% Convergence tolerance for eigs (tuning option). Smaller values give more
% accurate eigenvalues at the cost of more iterations.
tolerance = 1e-10;

% Maximum number of iterations for eigs (tuning option). Increase it if the
% solver fails to converge on harder domains.
max_iterations = 300;

% WARNING: The domain specified by the indicator function must be embedded
% into the rectangular bounding box specified above. This property is not
% checked in the code; it is your responsibility that this assertion is
% valid.

% Indicator function describing the domain
phi = @(x,y) indicator_rectangle(x, y, 0, pi/2, 0, pi/4);   % small rectangle [0, pi/2] x [0, pi/4]

% Assemble the masked DST-based Laplace operator and the grid/domain info
[Lop, info] = make_dst_laplace_op_batched(x_range, y_range, M, phi);

% Compute the first k eigenvalues (partial spectrum) with an iterative solver.
% 'largestreal' selects the eigenvalues closest to zero (least negative), i.e.
% the leading ones. The operator is symmetric, so IsFunctionSymmetric is set.
% Take the real part as a hedge against spurious imaginary parts.

D = eigs(Lop, info.dofs, k, 'largestreal', ...
    'Tolerance',           tolerance, ...
    'MaxIterations',       max_iterations, ...
    'SubspaceDimension',   subspace_dim, ...
    'IsFunctionSymmetric', true, ...
    'Display',             1);

eigs_numerical = sort(real(D), 'descend');
