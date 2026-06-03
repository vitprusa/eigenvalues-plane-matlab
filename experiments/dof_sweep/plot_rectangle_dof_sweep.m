function plot_rectangle_dof_sweep()
%PLOT_RECTANGLE_DOF_SWEEP Plot the rectangle DOF-sweep spectra from the CSVs.
%
%   Reads results/eigenvalues_dof_sweep/rectangle_<method>_*-eigenvalues.csv and
%   rectangle_analytic-eigenvalues.csv (produced by COMPUTE_RECTANGLE_DOF_SWEEP)
%   and saves rectangle_dof_sweep.png: eigenvalue index vs eigenvalue with all
%   method curves and the analytic ground truth. Both panels are clipped to the
%   physical band (the spurious high modes run off the top), and a lower-right
%   inset zooms to indices n <= 600.

    dof_sweep_dir   = fileparts(mfilename('fullpath'));
    experiments_dir = fileparts(dof_sweep_dir);
    addpath(dof_sweep_dir);     % load_dof_sweep
    addpath(experiments_dir);   % read_eigs_csv
    project_root = fileparts(experiments_dir);
    out_dir = fullfile(project_root, 'results', 'eigenvalues_dof_sweep');

    methods = {'dst', 'fd', 'fem', 'cheb'};
    results = load_dof_sweep(out_dir, 'rectangle', methods);

    % Largest DOF across all methods, used to truncate the analytic line.
    maxdof = 0;
    for mi = 1:numel(methods)
        maxdof = max(maxdof, max(cellfun(@(r) r.dofs, results.(methods{mi}))));
    end
    analytic = read_eigs_csv(fullfile(out_dir, 'rectangle_analytic-eigenvalues.csv'));
    analytic = analytic(1:min(maxdof, numel(analytic)));

    colors = struct('dst', [0 0.45 0.74], 'fd', [0.85 0.33 0.10], ...
                    'fem', [0.47 0.67 0.19], 'cheb', [0.49 0.18 0.56]);
    labels = struct('dst', 'DST', 'fd', 'FD', 'fem', 'FEM', 'cheb', 'Cheb');
    styles = {':', '-.', '--', '-'};

    % Render all text (labels, legend, tick labels) with the LaTeX
    % interpreter, i.e. in the standard LaTeX Computer Modern font.
    fig = figure('Visible', 'off', 'Position', [100 100 1000 720], ...
        'defaultAxesTickLabelInterpreter', 'latex', ...
        'defaultTextInterpreter',          'latex', ...
        'defaultLegendInterpreter',        'latex');

    main = axes(fig);
    draw_all(main, results, methods, analytic, colors, labels, styles, [], true);
    xlabel(main, 'eigenvalue index $n$');
    ylabel(main, '$\lambda_n$');
    title(main, 'Rectangle $[0, 2\pi] \times [0, \pi]$ --- eigenvalue DOF sweep');
    % One column per method (DST, FD, FEM, Cheb) plus a column for the analytic.
    legend(main, 'Location', 'northwest', 'NumColumns', numel(methods) + 1, 'FontSize', 7);
    grid(main, 'on'); box(main, 'on');
    % Clip y to the physical band; the spurious high modes (Cheb, FEM) run off
    % the top of the axes.
    ylim(main, [0, 1.5 * analytic(end)]);

    % Inset (lower-right, with a gap from the main axes): zoom to indices
    % n <= 600, with y clipped to the low (physical) eigenvalues.
    inset = axes(fig, 'Position', [0.58 0.15 0.304 0.304], 'Color', 'w');
    draw_all(inset, results, methods, analytic, colors, labels, styles, 600);
    grid(inset, 'on'); box(inset, 'on');
    ylim(inset, [0, 2 * analytic(min(600, numel(analytic)))]);
    title(inset, 'indices $n \leq 600$', 'FontSize', 8);
    set(inset, 'FontSize', 7);

    png = fullfile(out_dir, 'rectangle_dof_sweep.png');
    try
        exportgraphics(fig, png, 'Resolution', 150);
    catch
        print(fig, png, '-dpng', '-r150');
    end
    fprintf('Wrote plot %s\n', png);
end


function draw_all(ax, results, methods, analytic, colors, labels, styles, nmax, pad_legend)
    if nargin < 9
        pad_legend = false;
    end
    % Longest method, so the shorter ones can be padded to a full legend column.
    maxcount = max(cellfun(@(m) numel(results.(m)), methods));
    hold(ax, 'on');
    for mi = 1:numel(methods)
        m = methods{mi};
        nk = numel(results.(m));
        for k = 1:nk
            r = results.(m){k};
            if isempty(nmax)
                idx = 1:numel(r.evals);
            else
                idx = 1:min(nmax, numel(r.evals));
            end
            si = mod(k - 1, numel(styles)) + 1;               % cycle the 4 line styles
            lw = 1.2 + 0.9 * floor((k - 1) / numel(styles));  % thicker on each extra cycle
            plot(ax, idx, r.evals(idx), styles{si}, 'Color', colors.(m), 'LineWidth', lw, ...
                'DisplayName', sprintf('%s (dof = %d)', labels.(m), r.dofs));
        end
        % Pad this method's legend column to maxcount with invisible blank rows,
        % so the column-major legend keeps one column per method.
        if pad_legend
            for p = 1:(maxcount - nk)
                plot(ax, NaN, NaN, 'LineStyle', 'none', 'Marker', 'none', 'DisplayName', ' ');
            end
        end
    end
    % Analytic ground truth (black, thicker).
    if isempty(nmax)
        idx = 1:numel(analytic);
    else
        idx = 1:min(nmax, numel(analytic));
    end
    plot(ax, idx, analytic(idx), 'k-', 'LineWidth', 1.8, 'DisplayName', 'analytic (exact)');
    hold(ax, 'off');
    if ~isempty(nmax)
        xlim(ax, [0 nmax]);
    end
end
