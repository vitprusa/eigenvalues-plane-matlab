function [dst_laplace_op, info] = make_dst_laplace_op(x_range, y_range, M, indicator_function)
%MAKE_DST_LAPLACE_OP Factory for a masked DST-based discrete Laplace operator.
%
%   [dst_laplace_op, info] = make_dst_laplace_op(x_range, y_range, M, indicator_function)
%
%   Encapsulates the grid/domain set-up that is shared by all the
%   dst_laplace_grid_eigenvalues_*_test scripts and returns a single
%   operator handle, dst_laplace_op, that applies the discrete Laplacian
%   on the degrees of freedom (interior grid points selected by the
%   domain mask).
%
%   Inputs:
%     x_range            - 1-by-2 vector [x_min, x_max] giving the x-extent
%                          of the bounding box of the grid.
%     y_range            - 1-by-2 vector [y_min, y_max] giving the y-extent
%                          of the bounding box of the grid.
%     M                  - number of interior grid points along one full
%                          bounding-box slice;
%                          used to derive the (equal) grid spacing
%                          h = (x_max - x_min)/(M+1).
%     indicator_function - domain-defining function handle indicator_function(x, y) that is
%                          positive inside the domain. The mask keeps grid
%                          points with indicator_function(X, Y) > 0.
%
%   Outputs:
%     dst_laplace_op - function handle acting on a degrees-of-freedom
%                      column vector v (length info.dofs). It scatters v
%                      onto the grid using the mask, applies
%                      dst_laplace_grid, and gathers the result back into a
%                      degrees-of-freedom vector. The operator is symmetric
%                      (use 'IsFunctionSymmetric', true with eigs), even
%                      though an explicitly assembled matrix need not be.
%     info           - struct with the grid/domain data:
%                        .h            grid spacing
%                        .x_vec        grid coordinates in x
%                        .y_vec        grid coordinates in y
%                        .X, .Y        ndgrid coordinate arrays
%                        .global_mask  logical mask of interior grid points
%                        .dofs         number of degrees of freedom
%
%   Example:
%     indicator_function = @(x, y) indicator_rectangle(x, y, 0, pi, 0, pi);
%     [Lop, info] = make_dst_laplace_op([0 pi], [0 pi], 20, indicator_function);
%     L = zeros(info.dofs);
%     I = eye(info.dofs);
%     parfor i = 1:info.dofs
%         L(:, i) = Lop(I(:, i));
%     end
%     eigs_numerical = sort(real(eig(L)), 'descend');
%
%   See also DST_LAPLACE_GRID, VALS_VEC_TO_VALS_GRID, VALS_GRID_TO_VALS_VEC.

    % Grid spacing must be the same in both directions
    h = (x_range(2) - x_range(1)) / (M + 1);

    x_vec = x_range(1) : h : x_range(2);
    y_vec = y_range(1) : h : y_range(2);
    [X, Y] = ndgrid(x_vec, y_vec);

    % Create the domain mask from the defining function
    global_mask = indicator_function(X, Y) > 0;

    % Number of degrees of freedom (interior grid points after masking)
    dofs = sum(global_mask, 'all');

    % Laplace operator acting on a degrees-of-freedom values vector
    dst_laplace_op = @(v) vals_grid_to_vals_vec( ...
        dst_laplace_grid( ...
        vals_vec_to_vals_grid(v, global_mask), h), ...
        global_mask);

    info = struct( ...
        'h',           h, ...
        'x_vec',       x_vec, ...
        'y_vec',       y_vec, ...
        'X',           X, ...
        'Y',           Y, ...
        'global_mask', global_mask, ...
        'dofs',        dofs);

end
