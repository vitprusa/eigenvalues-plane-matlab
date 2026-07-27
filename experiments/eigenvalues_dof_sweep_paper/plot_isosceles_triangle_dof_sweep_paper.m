function plot_isosceles_triangle_dof_sweep_paper()
%PLOT_ISOSCELES_TRIANGLE_DOF_SWEEP_PAPER Paper variant of the triangle DOF-sweep plot.
%
%   Paper version of PLOT_ISOSCELES_TRIANGLE_DOF_SWEEP: same figure with no main
%   title, drawn in black and white. The method is encoded by line style (DST
%   dotted, FD dash-dot, FEM dashed) and the DOF level by a marker symbol
%   (circle, square, triangle, diamond, down-triangle for the successive DOFs),
%   with a few markers placed along each curve; the analytic ground truth is a
%   thick solid black line with no marker. Reads the existing CSVs from
%   results/eigenvalues_dof_sweep/ (produced by
%   COMPUTE_ISOSCELES_TRIANGLE_DOF_SWEEP) and writes
%   isosceles_triangle_dof_sweep.png into results_paper/eigenvalues_dof_sweep/.
%   The domain name is carried by the file name in place of the removed title.

    paper_dir       = fileparts(mfilename('fullpath'));
    experiments_dir = fileparts(paper_dir);
    orig_dir        = fullfile(experiments_dir, 'eigenvalues_dof_sweep');
    addpath(orig_dir);          % load_dof_sweep
    addpath(experiments_dir);   % read_eigs_csv
    project_root = fileparts(experiments_dir);
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

    % Black and white: method -> line style, DOF level -> marker symbol.
    mstyle  = struct('dst', ':', 'fd', '-.', 'fem', '--');
    markers = {'o', 's', '^', 'd', 'v'};
    labels  = struct('dst', 'DST', 'fd', 'FD', 'fem', 'FEM');

    % Render all text (labels, legend, tick labels) with the LaTeX
    % interpreter, i.e. in the standard LaTeX Computer Modern font.
    fig = figure('Visible', 'off', 'Position', [100 100 1000 720], ...
        'defaultAxesTickLabelInterpreter', 'latex', ...
        'defaultTextInterpreter',          'latex', ...
        'defaultLegendInterpreter',        'latex');

    main = axes(fig);
    draw_all(main, results, methods, analytic, mstyle, markers, labels, [], true);
    xlabel(main, 'eigenvalue index $k$');
    ylabel(main, '$\lambda_k$');
    % No title: the domain is identified by the output file name.
    % One column per method (DST, FD, FEM) plus a column for the analytic.
    legend(main, 'Location', 'northwest', 'NumColumns', numel(methods) + 1, 'FontSize', 7, 'Box', 'off');
    grid(main, 'on'); box(main, 'on');
    % Clip y to the physical band; the spurious high modes (FEM) run off the top.
    ylim(main, [0, 1.5 * analytic(end)]);

    % Inset (lower-right, with a gap from the main axes): zoom to indices
    % n <= 600, with y clipped to the low (physical) eigenvalues.
    inset = axes(fig, 'Position', [0.58 0.15 0.304 0.304], 'Color', 'w');
    draw_all(inset, results, methods, analytic, mstyle, markers, labels, 600);
    grid(inset, 'on'); box(inset, 'on');
    ylim(inset, [0, 2 * analytic(min(600, numel(analytic)))]);
    title(inset, 'indices $k \leq 600$', 'FontSize', 8);
    set(inset, 'FontSize', 7);

    png = fullfile(out_dir, 'isosceles_triangle_dof_sweep.png');
    try
        exportgraphics(fig, png, 'Resolution', 150);
    catch
        print(fig, png, '-dpng', '-r150');
    end
    fprintf('Wrote plot %s\n', png);
end


function draw_all(ax, results, methods, analytic, mstyle, markers, labels, nmax, pad_legend)
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
            % Method -> line style; DOF level -> marker symbol, drawn at a few
            % staggered points along the curve so the lines stay clean.
            mk = markers{mod(k - 1, numel(markers)) + 1};
            plot(ax, idx, r.evals(idx), mstyle.(m), 'Color', 'k', 'LineWidth', 1.1, ...
                'Marker', mk, 'MarkerSize', 5, 'MarkerEdgeColor', 'k', ...
                'MarkerFaceColor', 'none', 'MarkerIndices', marker_idx(numel(idx), k), ...
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
    % Analytic ground truth (solid black, thicker, no marker).
    if isempty(nmax)
        idx = 1:numel(analytic);
    else
        idx = 1:min(nmax, numel(analytic));
    end
    plot(ax, idx, analytic(idx), 'k-', 'LineWidth', 1.8, 'DisplayName', 'analytic (exact)');
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


function mi = marker_idx(n, k)
    % A handful of marker positions along a curve of length n, staggered by the
    % DOF index k so markers of overlapping curves do not all land together.
    nm = 3;
    if n <= 1
        mi = 1;
        return;
    end
    p = round(linspace(1, n, nm + 2));
    p = p(2:end-1);                       % drop the two endpoints
    shift = round((k - 1) / 5 * n / (nm + 1));
    mi = unique(min(n, max(1, p + shift)));
end
