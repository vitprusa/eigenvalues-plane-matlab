function [d2u_dx2, d2u_dy2, X, Y, mask] = dst_second_derivative_2d(phi, u_func, x_range, y_range, h)
% DST_SECOND_DERIVATIVE_2D  Computes the second derivative of a function u
% on a 2D domain defined implicitly by phi(x,y) <= 0, using the Discrete
% Sine Transform (DST) slice-by-slice.
%
% INPUTS:
%   phi     - function handle, implicit domain: phi(x,y) <= 0 defines interior
%             Example: @(x,y) x.^2 + y.^2 - 1  (unit disk)
%   u_func  - function handle u(x,y), the function to differentiate
%             Example: @(x,y) sin(pi*x) .* sin(pi*y)
%   x_range - [x_min, x_max], bounding box in x
%   y_range - [y_min, y_max], bounding box in y
%   h       - scalar grid spacing (same in both directions)
%
% OUTPUTS:
%   d2u_dx2 - matrix of d²u/dx² at grid points (NaN outside domain)
%   d2u_dy2 - matrix of d²u/dy² at grid points (NaN outside domain)
%   X, Y    - meshgrid coordinate matrices
%   mask    - logical matrix, true where point is inside the domain
%
% METHOD:
%   For each row (fixed y), the interior grid points form a 1D slice.
%   The DST-I is applied to that slice to obtain d²u/dx² spectrally.
%   The same is done column-wise for d²u/dy².
%
%   DST-I assumes homogeneous Dirichlet boundary conditions on the slice
%   endpoints. The second derivative in frequency space is:
%       d²u/dx² <=> -( k*pi/(N+1)/h )^2 * U_k
%   where k = 1,...,N and N is the number of interior points in the slice.
%
% EXAMPLE:
%   phi    = @(x,y) x.^2 + y.^2 - 0.9^2;   % disk of radius 0.9
%   u_func = @(x,y) (1 - x.^2 - y.^2);     % smooth, zero on boundary
%   [d2dx2, d2dy2, X, Y, mask] = dst_second_derivative_2d(phi, u_func, [-1,1], [-1,1], 0.05);
%   figure; surf(X, Y, d2dx2); title('d^2u/dx^2 via DST');

% -------------------------------------------------------------------------
% 1. Build the regular grid
% -------------------------------------------------------------------------
x_vec = x_range(1) : h : x_range(2);
y_vec = y_range(1) : h : y_range(2);
[X, Y] = meshgrid(x_vec, y_vec);          % size: (Ny x Nx)

% -------------------------------------------------------------------------
% 2. Evaluate u and build domain mask
% -------------------------------------------------------------------------
U    = u_func(X, Y);                       % function values on full grid
mask = phi(X, Y) <= 0;                     % true = inside domain

% Initialise output arrays with NaN (outside domain stays NaN)
d2u_dx2 = NaN(size(X));
d2u_dy2 = NaN(size(X));

% -------------------------------------------------------------------------
% 3. Row-wise DST: compute d²u/dx² for each row (fixed y_j)
% -------------------------------------------------------------------------
Ny = length(y_vec);

for j = 1 : Ny
    row_mask = mask(j, :);                 % 1 x Nx logical
    idx      = find(row_mask);             % column indices inside domain

    if numel(idx) < 2
        continue                           % skip rows with < 2 interior pts
    end

    N    = numel(idx);                     % number of points in this slice
    vals = U(j, idx);                      % u values on this slice

    % DST-I via MATLAB's dst (type-I sine transform)
    % dst_coeffs(k) corresponds to sin(k*pi*x / L_eff) basis,
    % where L_eff = (N+1)*h  (the effective length including ghost endpoints)
    dst_coeffs = dst1(vals);

    % Eigenvalues of the second-derivative operator under DST-I
    k_vec   = (1 : N)';                    % wave numbers
    L_eff   = (N + 1) * h;                % effective slice length
    lambda  = -( k_vec * pi / L_eff ).^2; % spectral multipliers

    % Multiply in frequency domain and invert
    d2_coeffs    = lambda' .* dst_coeffs;
    d2_vals      = idst1(d2_coeffs);

    d2u_dx2(j, idx) = d2_vals;
end

% -------------------------------------------------------------------------
% 4. Column-wise DST: compute d²u/dy² for each column (fixed x_i)
% -------------------------------------------------------------------------
Nx = length(x_vec);

for i = 1 : Nx
    col_mask = mask(:, i);                 % Ny x 1 logical
    idx      = find(col_mask);             % row indices inside domain

    if numel(idx) < 2
        continue
    end

    N    = numel(idx);
    vals = U(idx, i)';                     % row vector of u values

    dst_coeffs = dst1(vals);

    k_vec   = (1 : N)';
    L_eff   = (N + 1) * h;
    lambda  = -( k_vec * pi / L_eff ).^2;

    d2_coeffs      = lambda' .* dst_coeffs;
    d2_vals        = idst1(d2_coeffs);

    d2u_dy2(idx, i) = d2_vals';
end

end % ===== end main function ==============================================


% =========================================================================
% LOCAL HELPERS: DST-I and its inverse (identical up to normalisation)
% =========================================================================

function y = dst1(x)
% DST1  Type-I Discrete Sine Transform of row vector x.
%   Uses the FFT of an antisymmetric extension for efficiency.
%   Equivalent to MATLAB's dst(x, 'Type', 1) (Signal Processing Toolbox).
%   Falls back to a pure-FFT implementation so no toolbox is required.

N  = numel(x);
if N == 0, y = []; return; end

% Antisymmetric extension: [0, x, 0, -fliplr(x)]  length 2*(N+1)
ext = [0, x(:)', 0, -fliplr(x(:)')];
F   = fft(ext);

% Extract imaginary parts of elements 2 .. N+1 (scaled)
y = -imag( F(2 : N+1) );
end


function x = idst1(y)
% IDST1  Inverse of DST-I.  For DST-I the inverse is itself up to scaling.
%   IDST-I(y) = DST-I(y) * 2/(N+1)

N = numel(y);
if N == 0, x = []; return; end
x = dst1(y) * 2 / (N + 1);
end
