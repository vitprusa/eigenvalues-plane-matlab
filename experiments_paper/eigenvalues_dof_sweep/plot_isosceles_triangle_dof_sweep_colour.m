function plot_isosceles_triangle_dof_sweep_colour()
%PLOT_ISOSCELES_TRIANGLE_DOF_SWEEP_PAPER_COLOUR Colour paper variant of the triangle DOF sweep.
%
%   Colour counterpart of PLOT_ISOSCELES_TRIANGLE_DOF_SWEEP_PAPER: same figure
%   with no main title, but the method is encoded by colour (DST blue, FD orange,
%   FEM green) and the DOF level by line style (a per-cycle width bump keeps the
%   fifth DOF distinct from the first); lines are drawn thicker than in the
%   black-and-white variant. The analytic ground truth is a thick solid black
%   line. Reads the existing CSVs from results/eigenvalues_dof_sweep/ (produced
%   by COMPUTE_ISOSCELES_TRIANGLE_DOF_SWEEP) and writes
%   isosceles_triangle_dof_sweep_colour.eps into
%   results_paper/eigenvalues_dof_sweep/. The domain name is carried by the file
%   name in place of the removed title.

    paper_dir       = fileparts(mfilename('fullpath'));
    project_root    = fileparts(fileparts(paper_dir));
    experiments_dir = fullfile(project_root, 'experiments');
    orig_dir        = fullfile(experiments_dir, 'eigenvalues_dof_sweep');
    % Appended, not prepended: experiments/eigenvalues_dof_sweep holds
    % same-named plot functions, and prepending it would shadow the paper
    % ones for every call after the first in a batch run.
    addpath(orig_dir, '-end');          % load_dof_sweep
    addpath(experiments_dir, '-end');   % read_eigs_csv
    data_dir = fullfile(project_root, 'results', 'eigenvalues_dof_sweep');
    out_dir  = fullfile(project_root, 'results_paper', 'eigenvalues_dof_sweep');
    if ~exist(out_dir, 'dir')
        mkdir(out_dir);
    end

    methods = {'dst', 'fd', 'fem'};
    results = load_dof_sweep(data_dir, 'isosceles_triangle', methods);

    % Largest DOF across all methods, used to truncate the analytic line.
    maxdof = 0;
    for mi = 1:numel(methods)
        maxdof = max(maxdof, max(cellfun(@(r) r.dofs, results.(methods{mi}))));
    end
    analytic = read_eigs_csv(fullfile(data_dir, 'isosceles_triangle_analytic-eigenvalues.csv'));
    analytic = analytic(1:min(maxdof, numel(analytic)));

    % Colour: method -> colour, DOF level -> line style.
    colors = struct('dst', [0 0.45 0.74], 'fd', [0.85 0.33 0.10], ...
                    'fem', [0.47 0.67 0.19]);
    labels = struct('dst', 'DST', 'fd', 'FD', 'fem', 'FEM');
    styles = {':', '-.', '--', '-'};

    % Render all text (labels, legend, tick labels) with the LaTeX
    % interpreter, i.e. in the standard LaTeX Computer Modern font.
    fig = figure('Visible', 'off', 'Position', [100 100 1000 720], ...
        'defaultAxesTickLabelInterpreter', 'latex', ...
        'defaultTextInterpreter',          'latex', ...
        'defaultLegendInterpreter',        'latex');

    main = axes(fig);
    draw_all(main, results, methods, analytic, colors, labels, styles, [], true);
    xlabel(main, 'eigenvalue index $k$');
    ylabel(main, '$\lambda_k$');
    % No title: the domain is identified by the output file name.
    % One column per method (DST, FD, FEM) plus a column for the analytic.
    legend(main, 'Location', 'northwest', 'NumColumns', numel(methods) + 1, 'FontSize', 10.5, 'Box', 'off');
    grid(main, 'on'); box(main, 'on');
    % Clip y to the physical band; the spurious high modes (FEM) run off the top.
    ylim(main, [0, 1.5 * analytic(end)]);

    % Inset (lower-right, with a gap from the main axes): zoom to indices
    % n <= 600, with y clipped to the low (physical) eigenvalues.
    inset = axes(fig, 'Position', [0.58 0.15 0.304 0.304], 'Color', 'w');
    draw_all(inset, results, methods, analytic, colors, labels, styles, 600);
    grid(inset, 'on'); box(inset, 'on');
    ylim(inset, [0, 2 * analytic(min(600, numel(analytic)))]);
    title(inset, 'indices $k \leq 600$', 'FontSize', 8);
    set(inset, 'FontSize', 7);

    eps_file = fullfile(out_dir, 'isosceles_triangle_dof_sweep_colour.eps');
    try
        exportgraphics(fig, eps_file, 'ContentType', 'vector');
    catch
        print(fig, eps_file, '-depsc2', '-painters');
    end
    fprintf('Wrote plot %s\n', eps_file);
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
            lw = 2.0 + 1.0 * floor((k - 1) / numel(styles));  % thicker on each extra cycle
            plot(ax, idx, r.evals(idx), styles{si}, 'Color', colors.(m), 'LineWidth', lw, ...
                'DisplayName', sprintf('%s (DOF = %d)', labels.(m), r.dofs));
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
    plot(ax, idx, analytic(idx), 'k-', 'LineWidth', 2.6, 'DisplayName', 'analytic/exact');
    % Pad the analytic legend column to maxcount as well, so the column-major
    % legend has exactly maxcount rows and every method (plus analytic) keeps
    % its own column regardless of the method count.
    if pad_legend
        for p = 1:(maxcount - 1)
            plot(ax, NaN, NaN, 'LineStyle', 'none', 'Marker', 'none', 'DisplayName', ' ');
        end
    end
    hold(ax, 'off');
    if ~isempty(nmax)
        xlim(ax, [0 nmax]);
    end
end
