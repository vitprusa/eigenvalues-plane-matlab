% DST_LAPLACE_FULL_SPECTRUM_SMALL_RECTANGLE Full Dirichlet-Laplacian spectrum on a small rectangle.
%
%   Computes the complete set of eigenvalues of the discrete Laplace operator
%   with homogeneous Dirichlet boundary conditions on a small rectangle
%   [0, pi/2] x [0, pi/4] sitting inside the bounding box [0, pi] x [0, pi].
%   The domain is described by an indicator function embedded in a rectangular
%   bounding box and discretised on a uniform grid; the DST-based Laplace
%   matrix is assembled with make_dst_laplace_mat_batched and its full spectrum
%   is obtained with eig and sorted in descending order.
%
%   Generated from the dst_laplace_full_spectrum template.

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
M = 47;

% WARNING: The domain specified by the indicator function must be embedded
% into the rectangular bounding box specified above. This property is not
% checked in the code; it is your responsibility that this assertion is
% valid.

% Indicator function describing the domain
phi = @(x,y) indicator_rectangle(x, y, 0, pi/2, 0, pi/4);   % small rectangle [0, pi/2] x [0, pi/4]

% Assemble the masked DST-based Laplace matrix and the grid/domain info
[L, info] = make_dst_laplace_mat_batched(x_range, y_range, M, phi);

% Compute the eigenvalues (full spectrum)
% Take the real part as a hedge against spurious imaginary parts
[~, D] = eig(L);
eigs_numerical = sort(real(diag(D)), 'descend');
