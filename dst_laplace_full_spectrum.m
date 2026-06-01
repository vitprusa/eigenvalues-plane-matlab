% DST_LAPLACE_FULL_SPECTRUM Compute the full Dirichlet-Laplacian spectrum on a domain.
%
%   Script that computes the complete set of eigenvalues of the discrete
%   Laplace operator with homogeneous Dirichlet boundary conditions on a 2D
%   domain. The domain is described by an indicator function embedded in a
%   rectangular bounding box and discretised on a uniform grid; the DST-based
%   Laplace matrix is assembled with make_dst_laplace_mat_batched and its full
%   spectrum is obtained with eig and sorted in descending order.
%
%   Edit the bounding box (a, b, c, d), the resolution M, and the indicator
%   function phi below to change the domain and grid.

% Bounding box dimensions [a, b] x [c, d]
a = 0;
b = pi;
c = 0;
d = pi;

[x_range, y_range] = bounding_box(a, b, c, d);

% Number of DOF (interior grid points) in one slice, full bounding box
M = 20;

% WARNING: The domain specified by the indicator function must be embedded
% into the rectangular bounding box specified above. This property is not
% checked in the code; it is your responsibility that this assertion is
% valid.

% Indicator function describing the domain
phi = @(x,y) indicator_rectangle(x, y, a, b, c, d);   % rectangle

% Assemble the masked DST-based Laplace matrix and the grid/domain info
% [L, info] = make_dst_laplace_mat(x_range, y_range, M, phi);
[L, info] = make_dst_laplace_mat_batched(x_range, y_range, M, phi);

% Compute the eigenvalues (full spectrum)
% Take the real part as a hedge against spurious imaginary parts
[~, D] = eig(L);
eigs_numerical = sort(real(diag(D)), 'descend');

