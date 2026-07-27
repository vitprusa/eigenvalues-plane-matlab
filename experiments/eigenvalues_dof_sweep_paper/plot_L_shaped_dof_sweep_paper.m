function plot_L_shaped_dof_sweep_paper()
%PLOT_L_SHAPED_DOF_SWEEP_PAPER Paper variant of the L-shaped DOF-sweep plot.
%
%   Paper version of PLOT_L_SHAPED_DOF_SWEEP: same figure with no main title,
%   drawn in black and white. The method is encoded by line style (DST dotted,
%   FD dash-dot, FEM dashed) and the DOF level by a marker symbol (circle,
%   square, triangle, diamond for the successive DOFs), with a few markers
%   placed along each curve. Reads the existing CSVs from
%   results/eigenvalues_dof_sweep/ (produced by COMPUTE_L_SHAPED_DOF_SWEEP) and
%   writes L_shaped_dof_sweep.png into results_paper/eigenvalues_dof_sweep/. The
%   domain name is carried by the file name in place of the removed title.

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
    results = load_dof_sweep(data_dir, 'L_shaped', methods);

    % Black and white: method -> line style, DOF level -> marker symbol.
    mstyle  = struct('dst', ':', 'fd', '-.', 'fem', '--');
    markers = {'o', 's', '^', 'd', 'v'};
    labels  = struct('dst', 'DST', 'fd', 'FD', 'fem', 'FEM');

    % Render all text (labels, legend, tick labels) with the LaTeX
    % interpreter, i.e. in the standard LaTeX Computer Modern font.
    fig = figure('Visible', 'off', 'Position', [100 100 950 680], ...
        'defaultAxesTickLabelInterpreter', 'latex', ...
        'defaultTextInterpreter',          'latex', ...
        'defaultLegendInterpreter',        'latex');

    main = axes(fig);
    draw_curves(main, results, methods, mstyle, markers, labels, []);
    xlabel(main, 'eigenvalue index $k$');
    ylabel(main, '$\lambda_k$');
    % No title: the domain is identified by the output file name.
    legend(main, 'Location', 'southeast', 'NumColumns', 3, 'FontSize', 8, 'Box', 'off');
    grid(main, 'on'); box(main, 'on');

    inset = axes(fig, 'Position', [0.18 0.50 0.384 0.384], 'Color', 'w');
    draw_curves(inset, results, methods, mstyle, markers, labels, 1000);
    grid(inset, 'on'); box(inset, 'on');
    title(inset, 'indices $k \leq 1000$', 'FontSize', 8);
    set(inset, 'FontSize', 7);

    png = fullfile(out_dir, 'L_shaped_dof_sweep.png');
    try
        exportgraphics(fig, png, 'Resolution', 150);
    catch
        print(fig, png, '-dpng', '-r150');
    end
    fprintf('Wrote plot %s\n', png);
end


function draw_curves(ax, results, methods, mstyle, markers, labels, nmax)
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
                'DisplayName', sprintf('%s (dof = %d)', labels.(m), r.dofs));
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
