function laplace_vals_grid = dst_laplace_grid(vals_grid, h, num_workers)
%DST_LAPLACE_GRID Compute the Laplacian of a 2D grid using DST-based 1D second-derivative slices.
%
%   laplace_vals_grid = dst_laplace_grid(vals_grid, h)
%   laplace_vals_grid = dst_laplace_grid(vals_grid, h, num_workers)
%
%   Inputs:
%     vals_grid   - m-by-n numeric array of function values on a uniform grid.
%                   Each row corresponds to a fixed y and varying x, each column
%                   corresponds to a fixed x and varying y.
%     h           - scalar grid spacing (assumed equal in both x and y).
%     num_workers - (optional) number of workers for parfor. If omitted or Inf,
%                   uses the default parallel pool behavior.
%
%   Output:
%     laplace_vals_grid - m-by-n array containing the discrete Laplacian
%                         (d^2/dx^2 + d^2/dy^2) of vals_grid computed by
%                         applying a DST-based second-derivative routine
%                         to each column (x-derivative) and each row
%                         (y-derivative) and summing the results.
%
%   Notes:
%     - This function delegates 1D second-derivative computation to
%       dst_d2_slice, which is expected to accept a vector slice and the
%       spacing h and return the second derivative of that slice.
%     - The implementation uses parfor loops to compute derivatives along
%       columns and rows independently; set num_workers to control
%       parallelism.

    if nargin < 3
        num_workers = Inf;
    end

    vals_grid_D2XX = NaN(size(vals_grid));
    vals_grid_D2YY = NaN(size(vals_grid));

    [m, n] = size(vals_grid);

    % Compute the second derivative on column slices
    parfor (j = 1:n, num_workers)
        vals_grid_D2XX(:, j) = dst_d2_slice(vals_grid(:, j), h);
    end

    % Compute the second derivative on row slices
    parfor (j = 1:m, num_workers)
        vals_grid_D2YY(j, :) = dst_d2_slice(vals_grid(j, :), h);
    end

    laplace_vals_grid = vals_grid_D2XX + vals_grid_D2YY;

end
