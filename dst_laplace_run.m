clc;
clear

% Domain [a, b] x [c, d]
a = 0;
b = pi;
c = 0;
d = pi;

x_range = [a, b];
y_range = [c, d];

% Grid spacing must be the same in both directions
% M is the number of DOF (interior grid points) in one slice (before applying the mask)
M = 20;
h = (b-a)/(M+1);

x_vec = x_range(1) : h : x_range(2);
y_vec = y_range(1) : h : y_range(2);
[X, Y] = ndgrid(x_vec, y_vec); 

% Defining function for the domain
% phi = @(x,y) (x-pi/2).^2 + (y-pi/2).^2 - pi/2;   % unit disk
phi = @(x,y) indicator_rectangle(x, y, a, b, c, d);   % rectangle

% Create mask
global_mask = phi(X, Y) > 0; 

% Indexing, just for debugging
%seq = reshape(1:numel(X), size(X));

%masked_seq = NaN(size(X));
%masked_seq(global_mask) = seq(global_mask);
%masked_uvals_grid = masked_seq; 

% We need a function that vanishes on boundary
mm = 5;
nn = 1;
u = @(x, y) sin(mm*x).*sin(nn*y);
lapu = @(x, y) (-mm.^2 - nn.^2)*u(x, y)

u_grid = u(X, Y);
lapu_grid = lapu(X, Y);

masked_uvals_grid = NaN(size(X));
masked_uvals_grid(global_mask) = u_grid(global_mask);  

masked_lapu_grid = NaN(size(X));
masked_lapu_grid(global_mask) = lapu_grid(global_mask);  

[m, n] = size(X);

% Can I rewrite this as parfor loop?
valsgridD2XX = NaN(size(X));   
for j = 1:n
    slice_mask = global_mask(:, j);
    slice = masked_uvals_grid(slice_mask, j);
    if (isempty(slice))
        % Do nothing, the slice is empty.
    else
        % Take the second derivative on the slice.
        valsgridD2XX(slice_mask, j) = dst_d2_slice(slice, h);
    end
end

valsgridD2YY = NaN(size(X));   
for j = 1:m
    slice_mask = global_mask(j, :);
    slice = masked_uvals_grid(j, slice_mask);
    if (isempty(slice))
        % Do nothing, the slice is empty.
    else
        % Take the second derivative on the slice.
        valsgridD2YY(j, slice_mask) = dst_d2_slice(slice', h)';
    end
end

masked_lapu_grid;
laplace_vals_grid = valsgridD2XX + valsgridD2YY;

laplace_vals_grid_alt = dst_laplace_grid(masked_uvals_grid, global_mask, h);


residual = masked_lapu_grid - laplace_vals_grid
% Calculate the norm of the residual for convergence check
normResidual = norm(residual(global_mask));
disp(['Norm of the residual: ', num2str(normResidual)]);

% Plot function U on the masked domain
figure;
surf(X, Y, masked_uvals_grid, 'EdgeColor', 'none');
colormap(parula);
colorbar;
axis equal tight;
view(45, 30);
xlabel('x');
ylabel('y');
zlabel('U(x,y)');
title('Function U on the Masked Domain');
% Improve lighting for clarity
camlight headlight;
lighting gouraud;

% Visualize the residual on the masked domain
figure;
% Prepare Z for plotting: set outside points to NaN so they don't appear
Z = NaN(size(residual));
Z(global_mask) = residual(global_mask);

% Surface plot of residual
surf(X, Y, Z, 'EdgeColor', 'none');
colormap(jet);
colorbar;
axis equal tight;
view(45, 30);
xlabel('x');
ylabel('y');
zlabel('Residual (DU - \Delta U)');
title('Residual on the Masked Domain');
camlight headlight;
lighting gouraud;

% Also show a 2D color plot for clearer magnitude visualization
figure;
imagesc(x_vec, y_vec, Z.');
set(gca,'YDir','normal');
axis equal tight;
colormap(jet);
c = colorbar;
c.Label.String = 'Residual (DU - \Delta U)';
xlabel('x');
ylabel('y');
title('Residual (2D) on Masked Domain');

% Histogram of residual values inside the mask
figure;
histogram(residual(global_mask), 50);
xlabel('Residual value');
ylabel('Count');
title('Histogram of Residuals (inside mask)');

% Print max/min residual inside the mask
res_vals = residual(global_mask);
fprintf('Residual min: %g, max: %g, L2 norm: %g\n', min(res_vals), max(res_vals), norm(res_vals));