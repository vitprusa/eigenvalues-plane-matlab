function vals_grid = vals_vec_to_vals_grid(vals_vec, grid_mask)
    vals_grid = NaN(size(grid_mask));
    vals_grid(grid_mask) = vals_vec;
end