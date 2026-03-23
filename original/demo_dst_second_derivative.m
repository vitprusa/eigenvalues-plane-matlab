%% demo_dst_second_derivative.m
% Demonstrates dst_second_derivative_2d.m and validates the result
% against the analytical second derivative.
%
% Test case: unit disk,  u(x,y) = (1 - x^2 - y^2)
%   d²u/dx² = -2   (constant inside disk)
%   d²u/dy² = -2   (constant inside disk)

clear; clc; close all;

% ── Domain & function ────────────────────────────────────────────────────
r      = 0.9;
phi    = @(x,y)  x.^2 + y.^2 - r^2;          % disk, radius r
u_func = @(x,y)  (r^2 - x.^2 - y.^2);        % zero on boundary circle
d2u_analytical = -2;                           % both d²u/dx² and d²u/dy²

h = 0.02;                                      % grid spacing

% ── Call the function ─────────────────────────────────────────────────────
[d2dx2, d2dy2, X, Y, mask] = dst_second_derivative_2d( ...
        phi, u_func, [-1, 1], [-1, 1], h);

% ── Error analysis ────────────────────────────────────────────────────────
err_x = d2dx2(mask) - d2u_analytical;
err_y = d2dy2(mask) - d2u_analytical;

fprintf('=== DST Second Derivative – Verification ===\n');
fprintf('Grid spacing h = %.4f\n', h);
fprintf('Interior points: %d\n', sum(mask(:)));
fprintf('\nd²u/dx²:  max|error| = %.2e,  RMS error = %.2e\n', ...
        max(abs(err_x)), rms(err_x));
fprintf('d²u/dy²:  max|error| = %.2e,  RMS error = %.2e\n\n', ...
        max(abs(err_y)), rms(err_y));

% ── Plots ─────────────────────────────────────────────────────────────────
figure('Name','DST Second Derivative Demo','Position',[100 100 1100 380]);

subplot(1,3,1)
pcolor(X, Y, d2dx2); shading flat; axis equal tight; colorbar;
title('d²u/dx²  (DST)'); xlabel('x'); ylabel('y');

subplot(1,3,2)
pcolor(X, Y, d2dy2); shading flat; axis equal tight; colorbar;
title('d²u/dy²  (DST)'); xlabel('x'); ylabel('y');

subplot(1,3,3)
err_map = NaN(size(X));
err_map(mask) = abs(err_x);
pcolor(X, Y, err_map); shading flat; axis equal tight; colorbar;
title('|error|  d²u/dx²'); xlabel('x'); ylabel('y');

sgtitle(sprintf('Disk domain, h = %.3f', h));
