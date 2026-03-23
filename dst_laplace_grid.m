function laplace_vals_grid = dst_laplace_grid(vals_grid, grid_mask, h)

vals_grid_D2XX = NaN(size(grid_mask));
vals_grid_D2YY = NaN(size(grid_mask));

[m, n] = size(grid_mask);

% Compute the second derivative on row slices
parfor j = 1:n
    slice_mask = grid_mask(:, j);
    vals_slice = vals_grid(slice_mask, j);
    if (isempty(vals_slice))
        % Do nothing, the slice is empty.
    else
        % NaN vector
        d2_vals_slice = NaN(size(slice_mask))
        % Compute the derivative on the slice, NaN elements only
        d2_vals_slice(slice_mask) = dst_d2_slice(vals_slice, h);
        % Copy computed derivative values to the j-th column of the global
        % array
        vals_grid_D2XX(:, j) = d2_vals_slice;
        % We are in fact doing this, but the construction above allows us
        % to use parfor
        % vals_grid_D2XX(slice_mask, j) = dst_d2_slice(vals_slice, h);
    end
end

% Compute the second derivative on column slices
parfor j = 1:m
    slice_mask = grid_mask(j, :);
    vals_slice = vals_grid(j, slice_mask);
    if (isempty(vals_slice))
        % Do nothing, the slice is empty.
    else
        % NaN vector
        d2_vals_slice = NaN(size(slice_mask))
        % Compute the derivative on the slice, NaN elements only
        d2_vals_slice(slice_mask) = dst_d2_slice(vals_slice', h)';
        % Copy computed derivative values to the j-th column of the global
        % array
        vals_grid_D2YY(j, :) = d2_vals_slice;
        % We are in fact doing this, but the construction above allows us
        % to use parfor
        % vals_grid_D2YY(j, slice_mask) = dst_d2_slice(vals_slice', h)';
    end
end

laplace_vals_grid = vals_grid_D2XX + vals_grid_D2YY;

end