function laplace_vals_grid = dst_laplace_grid(vals_grid, h)

vals_grid_D2XX = NaN(size(vals_grid));
vals_grid_D2YY = NaN(size(vals_grid));

[m, n] = size(vals_grid);

% Compute the second derivative on column slices
parfor j = 1:n
    vals_grid_D2XX(:, j) = dst_d2_slice(vals_grid(:, j), h);
end

% Compute the second derivative on row slices
parfor j = 1:m
    vals_grid_D2YY(j, :) = dst_d2_slice(vals_grid(j, :), h);
end

laplace_vals_grid = vals_grid_D2XX + vals_grid_D2YY;

end
