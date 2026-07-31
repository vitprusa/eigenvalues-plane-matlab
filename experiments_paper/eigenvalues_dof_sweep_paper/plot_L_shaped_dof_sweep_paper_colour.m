function plot_L_shaped_dof_sweep_paper_colour()
%PLOT_L_SHAPED_DOF_SWEEP_PAPER_COLOUR Colour paper variant of the L-shaped DOF sweep.
%
%   Colour counterpart of PLOT_L_SHAPED_DOF_SWEEP_PAPER: same figure with no main
%   title, but the method is encoded by colour (DST blue, FD orange, FEM green)
%   and the DOF level by line style (a per-cycle width bump keeps the fifth DOF
%   distinct from the first); lines are drawn thicker than in the black-and-white
%   variant. Reads the existing CSVs from results/eigenvalues_dof_sweep/
%   (produced by COMPUTE_L_SHAPED_DOF_SWEEP) and writes
%   L_shaped_dof_sweep_colour.eps into results_paper/eigenvalues_dof_sweep/. The
%   domain name is carried by the file name in place of the removed title.

    paper_dir       = fileparts(mfilename('fullpath'));
    project_root    = fileparts(fileparts(paper_dir));
    experiments_dir = fullfile(project_root, 'experiments');
    orig_dir        = fullfile(experiments_dir, 'eigenvalues_dof_sweep');
    addpath(orig_dir);          % load_dof_sweep
    addpath(experiments_dir);   % read_eigs_csv
    data_dir = fullfile(project_root, 'results', 'eigenvalues_dof_sweep');
    out_dir  = fullfile(project_root, 'results_paper', 'eigenvalues_dof_sweep');
    if ~exist(out_dir, 'dir')
        mkdir(out_dir);
    end

    methods = {'dst', 'fd', 'fem'};
    results = load_dof_sweep(data_dir, 'L_shaped', methods);

    % Colour: method -> colour, DOF level -> line style.
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
    xlabel(main, 'eigenvalue index $k$');
    ylabel(main, '$\lambda_k$');
    % No title: the domain is identified by the output file name.
    legend(main, 'Location', 'southeast', 'NumColumns', 3, 'FontSize', 8, 'Box', 'off');
    grid(main, 'on'); box(main, 'on');

    inset = axes(fig, 'Position', [0.18 0.50 0.384 0.384], 'Color', 'w');
    draw_curves(inset, results, methods, colors, labels, styles, 1000);
    grid(inset, 'on'); box(inset, 'on');
    title(inset, 'indices $k \leq 1000$', 'FontSize', 8);
    set(inset, 'FontSize', 7);

    eps_file = fullfile(out_dir, 'L_shaped_dof_sweep_colour.eps');
    try
        exportgraphics(fig, eps_file, 'ContentType', 'vector');
    catch
        print(fig, eps_file, '-depsc2', '-painters');
    end
    fprintf('Wrote plot %s\n', eps_file);
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
            lw = 2.0 + 1.0 * floor((k - 1) / numel(styles));
            plot(ax, idx, r.evals(idx), styles{si}, 'Color', colors.(m), 'LineWidth', lw, ...
                'DisplayName', sprintf('%s (DOF = %d)', labels.(m), r.dofs));
        end
    end
    hold(ax, 'off');
    if ~isempty(nmax)
        xlim(ax, [0 nmax]);
    end
end
