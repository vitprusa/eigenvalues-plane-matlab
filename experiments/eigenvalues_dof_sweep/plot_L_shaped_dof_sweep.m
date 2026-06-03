function plot_L_shaped_dof_sweep()
%PLOT_L_SHAPED_DOF_SWEEP Plot the L-shaped DOF-sweep spectra from the CSVs.
%
%   Reads results/eigenvalues_dof_sweep/L_shaped_<method>_*-eigenvalues.csv
%   (produced by COMPUTE_L_SHAPED_DOF_SWEEP) and saves L_shaped_dof_sweep.png:
%   eigenvalue index vs eigenvalue with all curves (full spectra), plus an
%   upper-left inset zooming to indices n <= 1000.

    dof_sweep_dir   = fileparts(mfilename('fullpath'));
    experiments_dir = fileparts(dof_sweep_dir);
    addpath(dof_sweep_dir);     % load_dof_sweep
    addpath(experiments_dir);   % read_eigs_csv
    project_root = fileparts(experiments_dir);
    out_dir = fullfile(project_root, 'results', 'eigenvalues_dof_sweep');

    methods = {'dst', 'fd', 'fem'};
    results = load_dof_sweep(out_dir, 'L_shaped', methods);

    colors = struct('dst', [0 0.45 0.74], 'fd', [0.85 0.33 0.10], 'fem', [0.47 0.67 0.19]);
    labels = struct('dst', 'DST', 'fd', 'FD', 'fem', 'FEM');
    styles = {':', '-.', '--', '-'};

    % Render all text (labels, legend, tick labels) with the LaTeX
    % interpreter, i.e. in the standard LaTeX Computer Modern font.
    fig = figure('Visible', 'off', 'Position', [100 100 950 680], ...
        'defaultAxesTickLabelInterpreter', 'latex', ...
        'defaultTextInterpreter',          'latex', ...
        'defaultLegendInterpreter',        'latex');

    main = axes(fig);
    draw_curves(main, results, methods, colors, labels, styles, []);
    xlabel(main, 'eigenvalue index $n$');
    ylabel(main, '$\lambda_n$');
    title(main, 'L-shaped domain --- eigenvalue DOF sweep (DST, FD, FEM)');
    legend(main, 'Location', 'southeast', 'NumColumns', 3, 'FontSize', 8, 'Box', 'off');
    grid(main, 'on'); box(main, 'on');

    inset = axes(fig, 'Position', [0.18 0.50 0.384 0.384], 'Color', 'w');
    draw_curves(inset, results, methods, colors, labels, styles, 1000);
    grid(inset, 'on'); box(inset, 'on');
    title(inset, 'indices $n \leq 1000$', 'FontSize', 8);
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
            si = mod(k - 1, numel(styles)) + 1;
            lw = 1.2 + 0.9 * floor((k - 1) / numel(styles));
            plot(ax, idx, r.evals(idx), styles{si}, 'Color', colors.(m), 'LineWidth', lw, ...
                'DisplayName', sprintf('%s (dof = %d)', labels.(m), r.dofs));
        end
    end
    hold(ax, 'off');
    if ~isempty(nmax)
        xlim(ax, [0 nmax]);
    end
end
