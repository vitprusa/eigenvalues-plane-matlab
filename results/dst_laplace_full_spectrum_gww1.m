% DST_LAPLACE_FULL_SPECTRUM_GWW1 Full Dirichlet-Laplacian spectrum on the GWW1 isospectral drum.
%
%   Computes the complete set of eigenvalues of the discrete Laplace operator
%   with homogeneous Dirichlet boundary conditions on the first Gordon-Webb-
%   Wolpert (GWW1) isospectral drum. The domain is described by an indicator
%   function embedded in a rectangular bounding box and discretised on a
%   uniform grid; the DST-based Laplace matrix is assembled with
%   make_dst_laplace_mat_batched and its full spectrum is obtained with eig and
%   sorted in descending order.
%
%   GWW1 and GWW2 are isospectral, so this spectrum should match that of
%   dst_laplace_*_gww2.
%
%   Generated from the dst_laplace_full_spectrum template.

clc;
clear

% This script lives in results/; put the project sources on the path.
run(fullfile(fileparts(fileparts(mfilename('fullpath'))), 'startup.m'));

% Bounding box dimensions [a, b] x [c, d]
a = -3;
b =  3;
c = -3;
d =  3;

[x_range, y_range] = bounding_box(a, b, c, d);

% Number of DOF (interior grid points) in one slice, full bounding box
M = 50;

% WARNING: The domain specified by the indicator function must be embedded
% into the rectangular bounding box specified above. This property is not
% checked in the code; it is your responsibility that this assertion is
% valid.

% Indicator function describing the domain
phi = @(x,y) indicator_gww1(x, y, -3, -1, 1, 3, -1, 1);   % GWW1 isospectral drum

% Assemble the masked DST-based Laplace matrix and the grid/domain info
[L, info] = make_dst_laplace_mat_batched(x_range, y_range, M, phi);

% Compute the eigenvalues (full spectrum)
% Take the real part as a hedge against spurious imaginary parts
[~, D] = eig(L);
eigs_numerical = sort(real(diag(D)), 'descend');
