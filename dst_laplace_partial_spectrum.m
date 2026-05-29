% DST_LAPLACE_PARTIAL_SPECTRUM Compute the first few Dirichlet-Laplacian eigenvalues on a domain.
%
%   Script that computes only the first few eigenvalues of the discrete
%   Laplace operator with homogeneous Dirichlet boundary conditions on a 2D
%   domain. The domain is described by an indicator function embedded in a
%   rectangular bounding box and discretised on a uniform grid; the DST-based
%   Laplace operator is built with make_dst_laplace_op_batched and the leading
%   eigenvalues are obtained with the iterative solver eigs, sorted in
%   descending order.
%
%   This is the partial-spectrum counterpart of dst_laplace_full_spectrum:
%   use this when only the leading eigenvalues are needed and forming/solving
%   the full dense matrix would be wasteful.
%
%   Edit the bounding box (a, b, c, d), the resolution M, the number of
%   eigenvalues k, and the indicator function phi below to change the
%   problem.

clc;
clear

% Bounding box dimensions [a, b] x [c, d]
a = 0;
b = pi;
c = 0;
d = pi;

[x_range, y_range] = bounding_box(a, b, c, d);

% Number of DOF (interior grid points) in one slice, full bounding box
M = 200;

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
phi = @(x,y) indicator_rectangle(x, y, a, b, c, d);   % rectangle

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
