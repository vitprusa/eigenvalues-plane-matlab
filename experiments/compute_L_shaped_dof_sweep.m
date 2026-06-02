function compute_L_shaped_dof_sweep()
%COMPUTE_L_SHAPED_DOF_SWEEP DOF sweep of the L-shaped spectrum (DST, FD, FEM).
%
%   Computes the full Dirichlet-Laplacian spectrum on the L-shaped domain by
%   three methods -- DST, finite differences, and finite elements (all dense
%   eig) -- at four resolutions (degree-of-freedom settings) each, stores the
%   spectra as CSVs in results/eigenvalues_dof_sweep/, and saves a plot of the
%   eigenvalue index vs the eigenvalue with all 12 curves in one figure, plus an
%   inset zooming into the low end (indices n <= 1000).
%
%   Requires the PDE Toolbox (for the FEM part). Run headless with
%   experiments/run_L_shaped_dof_sweep.sh.

    experiments_dir = fileparts(mfilename('fullpath'));
    project_root    = fileparts(experiments_dir);
    run(fullfile(project_root, 'startup.m'));

    out_dir = fullfile(project_root, 'results', 'eigenvalues_dof_sweep');
    if ~exist(out_dir, 'dir')
        mkdir(out_dir);
    end

    % L-shaped domain: [-1,1]^2 minus the upper-right unit square [0,1]^2.
    box = [-1 1 -1 1];
    phi = @(x, y) indicator_L_shaped(x, y, -1, 0, 1, -1, 0, 1);   % DST / FD
    [x_range, y_range] = bounding_box(box(1), box(2), box(3), box(4));
    % FEM decsg geometry for the same L-shape (R1 - R2).
    gd = [[3; 4; -1; 1; 1; -1; -1; -1; 1; 1], [3; 4; 0; 1; 1; 0; 0; 0; 1; 1]];
    ns = char('R1', 'R2')';
    sf = 'R1-R2';

    % Four DOF settings per method.
    M_grid   = [15 25 35 50];          % DST and FD grid resolutions
    Hmax_fem = [0.20 0.13 0.09 0.06];  % FEM target mesh sizes
    nset = 4;

    results = struct('dst', {cell(1, nset)}, 'fd', {cell(1, nset)}, 'fem', {cell(1, nset)});

    % --- DST (full spectrum via dense eig) ---
    for k = 1:nset
        M = M_grid(k);
        [L, info] = make_dst_laplace_mat_batched(x_range, y_range, M, phi);
        evals = sort(-real(eig(L)), 'ascend');
        results.dst{k} = pack(evals, info.dofs, sprintf('M = %d', M));
    end

    % --- FD (5-point stencil, dense eig) ---
    for k = 1:nset
        M = M_grid(k);
        [evals, info] = fd_laplace_spectrum(struct('box', box, 'phi', phi, 'M', M));
        results.fd{k} = pack(evals, info.dofs, sprintf('M = %d', M));
    end

    % --- FEM (assembled K, M + dense eig) ---
    for k = 1:nset
        h = Hmax_fem(k);
        entry = struct('gd', gd, 'ns', ns, 'sf', sf, 'Hmax_eig', h, 'Hmax_solvepdeeig', h);
        [evals, info] = fem_laplace_spectrum(entry, "eig");
        results.fem{k} = pack(evals, info.dofs, sprintf('Hmax = %g', h));
    end

    % --- Store the spectra as CSVs ---
    methods = {'dst', 'fd', 'fem'};
    for mi = 1:numel(methods)
        m = methods{mi};
        for k = 1:nset
            r = results.(m){k};
            csv_file = fullfile(out_dir, sprintf('L_shaped_%s_%d-eigenvalues.csv', m, k));
            write_sweep_csv(csv_file, r, m);
        end
    end

    % --- Plot all 12 curves ---
    plot_sweep(results, methods, out_dir);
end


function r = pack(evals, dofs, res)
    r = struct('evals', evals(:), 'dofs', dofs, 'res', res);
end


function write_sweep_csv(csv_file, r, method)
    fid = fopen(csv_file, 'w');
    if fid == -1
        error('compute_L_shaped_dof_sweep:cannotOpen', 'Could not open %s.', csv_file);
    end
    fprintf(fid, '# Domain: L_shaped (%s, dense eig full spectrum)\n', upper(method));
    fprintf(fid, '# Computed %s\n', datestr(now, 'yyyy-mm-dd HH:MM:SS'));
    fprintf(fid, '# Resolution %s, dofs = %d\n', r.res, r.dofs);
    fclose(fid);
    n        = (1:numel(r.evals))';
    lambda_n = r.evals;
    writetable(table(n, lambda_n), csv_file, 'WriteMode', 'append', 'WriteVariableNames', true);
end


function plot_sweep(results, methods, out_dir)
    % Plot the full computed spectrum of every (method, resolution): all
    % eigenvalues, so each curve runs out to its own DOF count. An inset zooms
    % into the low end (indices n <= 1000) to show the behaviour near zero.
    colors = struct('dst', [0 0.45 0.74], 'fd', [0.85 0.33 0.10], 'fem', [0.47 0.67 0.19]);
    labels = struct('dst', 'DST', 'fd', 'FD', 'fem', 'FEM');
    styles = {':', '-.', '--', '-'};     % coarse -> fine DOF

    fig = figure('Visible', 'off', 'Position', [100 100 950 680]);

    main = axes(fig);
    draw_curves(main, results, methods, colors, labels, styles, []);
    xlabel(main, 'eigenvalue index n');
    ylabel(main, '\lambda_n');
    title(main, 'L-shaped domain — eigenvalue DOF sweep (DST, FD, FEM)');
    legend(main, 'Location', 'southeast', 'NumColumns', 3, 'FontSize', 8);
    grid(main, 'on'); box(main, 'on');

    % Inset (upper-left, clear of the main axes): zoom to indices n <= 1000.
    inset = axes(fig, 'Position', [0.18 0.50 0.384 0.384], 'Color', 'w');
    draw_curves(inset, results, methods, colors, labels, styles, 1000);
    grid(inset, 'on'); box(inset, 'on');
    title(inset, 'indices n \leq 1000', 'FontSize', 8);
    set(inset, 'FontSize', 7);

    png = fullfile(out_dir, 'L_shaped_dof_sweep.png');
    try
        exportgraphics(fig, png, 'Resolution', 150);
    catch
        print(fig, png, '-dpng', '-r150');
    end
    fprintf('Wrote plot %s\n', png);
end


function draw_curves(ax, results, methods, colors, labels, styles, nmax)
    % Plot every (method, resolution) curve on axes ax. If nmax is non-empty,
    % restrict to the first nmax eigenvalue indices and set that as the x-limit.
    hold(ax, 'on');
    for mi = 1:numel(methods)
        m = methods{mi};
        for k = 1:numel(results.(m))
            r = results.(m){k};
            if isempty(nmax)
                idx = 1:numel(r.evals);
            else
                idx = 1:min(nmax, numel(r.evals));
            end
            plot(ax, idx, r.evals(idx), styles{k}, 'Color', colors.(m), 'LineWidth', 1.2, ...
                'DisplayName', sprintf('%s (dof = %d)', labels.(m), r.dofs));
        end
    end
    hold(ax, 'off');
    if ~isempty(nmax)
        xlim(ax, [0 nmax]);
    end
end
