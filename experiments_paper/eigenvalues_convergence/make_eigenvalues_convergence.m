function make_eigenvalues_convergence(name)
%MAKE_EIGENVALUES_CONVERGENCE Ground-state convergence against DOF (paper).
%
%   For each domain it computes the first eigenvalue lambda_1 of the Dirichlet
%   Laplacian with DST, finite differences and finite elements over a sweep of
%   resolutions, and plots lambda_1 against the degrees of freedom actually used
%   by each discretisation. Unlike the DOF sweeps of
%   experiments_paper/eigenvalues_dof_sweep_paper -- which show the whole
%   spectrum at a few fixed resolutions -- this one follows a single eigenvalue
%   across many resolutions, so the horizontal axis is the DOF count itself.
%
%   The DST/FD grid resolutions M are chosen so that the DOF counts land near
%   200, 700, ..., 3000 while keeping the domain edges on grid lines; the FEM
%   mesh sizes Hmax target the same DOF counts. The three methods do not reach
%   identical DOF counts, which is why the measured count, not the resolution
%   parameter, is what gets plotted.
%
%   Output, into results_paper/eigenvalues_convergence/:
%     - <domain>_convergence.csv   one row per (method, resolution) with
%       the resolution parameter, the dof count, lambda_1 and the timing, and
%     - <domain>_convergence.png   the error |lambda_1 - lambda_1_ref|
%       against DOF on log-log axes, one curve per method, so that an algebraic
%       convergence rate shows up as a straight line of that slope. No title:
%       the domain is carried by the file name, as elsewhere in
%       experiments_paper.
%
%   An existing CSV is reused rather than recomputed, so the figure can be
%   restyled without paying for the spectra again; delete the CSV to force a
%   recompute.
%
%   Requires the PDE Toolbox (FEM).

    here         = fileparts(mfilename('fullpath'));
    project_root = fileparts(fileparts(here));
    run(fullfile(project_root, 'startup.m'));
    addpath(fullfile(project_root, 'experiments'));   % domain_catalog_dst / _fem

    out_dir = fullfile(project_root, 'results_paper', 'eigenvalues_convergence');
    if ~exist(out_dir, 'dir')
        mkdir(out_dir);
    end

    cfgs = domain_configs(project_root);
    if nargin >= 1 && ~isempty(name)
        cfgs = cfgs(strcmp({cfgs.name}, name));
        if isempty(cfgs)
            error('eigenvalues_convergence:unknownDomain', ...
                'Unknown domain "%s" (known: L_shaped).', name);
        end
    end
    for i = 1:numel(cfgs)
        cfg = cfgs(i);
        fprintf('=== %s ===\n', cfg.name);
        csv = fullfile(out_dir, sprintf('%s_convergence.csv', cfg.name));
        if exist(csv, 'file')
            runs = read_runs_csv(csv);
            fprintf('  (reusing %s)\n', csv);
        else
            runs = compute_runs(cfg);
            write_runs_csv(csv, cfg, runs);
            fprintf('  Wrote %s\n', csv);
        end
        png = fullfile(out_dir, sprintf('%s_convergence.png', cfg.name));
        plot_runs(png, cfg, runs);
        fprintf('  Wrote %s\n', png);
    end
end


function cfgs = domain_configs(project_root)
%DOMAIN_CONFIGS Geometry, resolution sweeps and reference lambda_1 per domain.
    cfgs = struct('name', {}, 'pretty', {}, 'box', {}, 'phi', {}, ...
                  'gd', {}, 'ns', {}, 'sf', {}, 'M_grid', {}, 'Hmax_fem', {}, ...
                  'lambda1_ref', {}, 'ref_label', {});

    % --- L-shaped domain ----------------------------------------------------
    % Box [-1,1]^2, so h = 2/(M+1) and M+1 must be even for the re-entrant
    % corner at the origin to fall on a grid line -- alignment matters a lot
    % here, the corner carrying the leading singularity. The M values sweep the
    % DOF counts from about 100 to about 10000: steps of roughly 500 up to 5000,
    % then coarser ones, the dense eig cost growing as dofs^3. The Hmax values
    % target the same counts through the empirical dofs ~ 12.8/Hmax^2.
    % Reference: MPS value of Betcke & Trefethen, results/eigenvalues/mps/.
    cfgs(end+1) = struct( ...
        'name', 'L_shaped', 'pretty', 'L-shaped', ...
        'box', [-1 1 -1 1], ...
        'phi', @(x, y) indicator_L_shaped(x, y, -1, 0, 1, -1, 0, 1), ...
        'gd', [[3; 4; -1; 1; 1; -1; -1; -1; 1; 1], [3; 4; 0; 1; 1; 0; 0; 0; 1; 1]], ...
        'ns', char('R1', 'R2')', 'sf', 'R1-R2', ...
        'M_grid',   [13 17 31 41 49 55 63 69 75 81 93 105 115], ...
        'Hmax_fem', [0.33 0.25 0.135 0.10 0.085 0.0755 0.066 0.060 0.055 0.051 ...
                     0.0446 0.0395 0.036], ...
        'lambda1_ref', 9.6397238440219, 'ref_label', 'MPS');
end


function runs = compute_runs(cfg)
%COMPUTE_RUNS One timed lambda_1 per method and resolution.
    [x_range, y_range] = bounding_box(cfg.box(1), cfg.box(2), cfg.box(3), cfg.box(4));

    runs = struct('method', {}, 'resolution', {}, 'dofs', {}, 'lambda1', {}, 'time', {});

    for k = 1:numel(cfg.M_grid)
        M = cfg.M_grid(k);

        t = tic;
        [L, info] = make_dst_laplace_mat_batched(x_range, y_range, M, cfg.phi);
        evals = sort(-real(eig(L)), 'ascend');
        runs(end+1) = mk('DST', sprintf('M = %d', M), info.dofs, evals(1), toc(t)); %#ok<AGROW>
        report(runs(end));

        t = tic;
        [evals, info] = fd_laplace_spectrum(struct('box', cfg.box, 'phi', cfg.phi, 'M', M));
        evals = sort(real(evals(:)), 'ascend');
        runs(end+1) = mk('FD', sprintf('M = %d', M), info.dofs, evals(1), toc(t)); %#ok<AGROW>
        report(runs(end));
    end

    for k = 1:numel(cfg.Hmax_fem)
        h = cfg.Hmax_fem(k);
        entry = struct('gd', cfg.gd, 'ns', cfg.ns, 'sf', cfg.sf, ...
                       'Hmax_eig', h, 'Hmax_solvepdeeig', h);
        t = tic;
        [evals, info] = fem_laplace_spectrum(entry, "eig");
        evals = sort(real(evals(:)), 'ascend');
        runs(end+1) = mk('FEM', sprintf('Hmax = %g', h), info.dofs, evals(1), toc(t)); %#ok<AGROW>
        report(runs(end));
    end
end


function s = mk(method, resolution, dofs, lambda1, tsec)
    s = struct('method', method, 'resolution', resolution, 'dofs', dofs, ...
               'lambda1', lambda1, 'time', tsec);
end


function report(r)
    fprintf('  %-4s %-12s dofs = %-5d  lambda_1 = %.8f  time = %6.2f s\n', ...
            r.method, r.resolution, r.dofs, r.lambda1, r.time);
end


function write_runs_csv(csv, cfg, runs)
%WRITE_RUNS_CSV One row per run: method, resolution, dofs, lambda_1, time.
    fid = fopen(csv, 'w');
    if fid == -1
        error('eigenvalues_convergence:csv', 'Could not open %s for writing.', csv);
    end
    closer = onCleanup(@() fclose(fid));
    fprintf(fid, '# Domain: %s (ground-state convergence against DOF, dense eig)\n', cfg.name);
    fprintf(fid, '# Reference lambda_1 (%s): %.13f\n', cfg.ref_label, cfg.lambda1_ref);
    fprintf(fid, 'method,resolution,dofs,lambda_1,time_s\n');
    for i = 1:numel(runs)
        fprintf(fid, '%s,%s,%d,%.12g,%.4f\n', runs(i).method, runs(i).resolution, ...
                runs(i).dofs, runs(i).lambda1, runs(i).time);
    end
end


function runs = read_runs_csv(csv)
%READ_RUNS_CSV Read back a run table written by WRITE_RUNS_CSV.
    fid = fopen(csv, 'r');
    if fid == -1
        error('eigenvalues_convergence:csvread', 'Could not open %s.', csv);
    end
    closer = onCleanup(@() fclose(fid));
    runs = struct('method', {}, 'resolution', {}, 'dofs', {}, 'lambda1', {}, 'time', {});
    while true
        line = fgetl(fid);
        if ~ischar(line); break; end
        s = strtrim(line);
        if isempty(s) || s(1) == '#' || strncmpi(s, 'method,', 7); continue; end
        p = strsplit(s, ',');
        if numel(p) < 5; continue; end
        runs(end+1) = mk(p{1}, p{2}, str2double(p{3}), ...
                         str2double(p{4}), str2double(p{5})); %#ok<AGROW>
    end
end


function plot_runs(png, cfg, runs)
%PLOT_RUNS Error in lambda_1 against DOF on log-log axes, one curve per method.
    methods = {'DST', 'FD', 'FEM'};
    colors  = {[0 0.45 0.74], [0.85 0.33 0.10], [0.47 0.67 0.19]};
    markers = {'o', 's', '^'};

    % Render all text with the LaTeX interpreter, as in the other paper figures.
    fig = figure('Visible', 'off', 'Position', [100 100 950 680], ...
        'defaultAxesTickLabelInterpreter', 'latex', ...
        'defaultTextInterpreter',          'latex', ...
        'defaultLegendInterpreter',        'latex');
    ax = axes(fig);
    hold(ax, 'on');

    % Each triangle sits above its own curve, staggered left to right so that no
    % two stack up and none reaches the curve above it. The DST-FD gap closes as
    % the grids refine, so FD takes the coarse left end where that gap is still
    % wide; DST then moves right of it, into the free space above every curve,
    % and FEM takes the fine right end where it has pulled away from FD.
    spans   = {[0.32 0.55], [0.08 0.28], [0.68 0.91]};
    offsets = [1.13, 1.25, 1.13];

    lo = inf; hi = 0;
    for mi = 1:numel(methods)
        sel = strcmp({runs.method}, methods{mi});
        d = [runs(sel).dofs];
        v = [runs(sel).lambda1];
        [d, ord] = sort(d);
        err = abs(v(ord) - cfg.lambda1_ref);
        % Least-squares algebraic rate: err ~ dofs^p, annotated by the triangle.
        p = polyfit(log(d), log(err), 1);
        loglog(ax, d, err, ['-' markers{mi}], 'Color', colors{mi}, 'LineWidth', 2.0, ...
            'MarkerSize', 7, 'MarkerFaceColor', 'w', ...
            'DisplayName', sprintf('\\texttt{%s}', methods{mi}));
        yhi = slope_triangle(ax, d, p, colors{mi}, spans{mi}, offsets(mi));
        lo = min([lo, err]);
        hi = max([hi, err, yhi]);
    end

    hold(ax, 'off');
    set(ax, 'XScale', 'log', 'YScale', 'log');
    ylim(ax, [lo * 0.72, hi * 1.35]);
    xlabel(ax, 'degrees of freedom');
    ylabel(ax, sprintf('$|\\lambda_1 - \\lambda_1^{\\mathrm{%s}}|$', cfg.ref_label));
    % No title: the domain is identified by the output file name.
    legend(ax, 'Location', 'southwest', 'FontSize', 10, 'Box', 'off');
    grid(ax, 'on'); box(ax, 'on');
    set(ax, 'FontSize', 12);

    try
        exportgraphics(fig, png, 'Resolution', 150);
    catch
        print(fig, png, '-dpng', '-r150');
    end
    close(fig);
end


function y_hi = slope_triangle(ax, d, p, color, span, offset)
%SLOPE_TRIANGLE Reference-slope triangle above one curve, labelled with the rate.
%
%   Right triangle whose hypotenuse has the fitted slope p(1), drawn just above
%   the curve it belongs to: horizontal leg on top, vertical leg on the right.
%   SPAN gives the fraction of the log-DOF range it covers and OFFSET how far
%   above the fitted line it sits -- the caller staggers the three triangles so
%   that each keeps to the gap above its own curve without running into the
%   neighbouring one. Returns the highest ordinate drawn so the caller can
%   leave room. Kept out of the legend.

    lg = log(d([1 end]));
    x1 = exp(lg(1) + span(1) * diff(lg));
    x2 = exp(lg(1) + span(2) * diff(lg));
    y1 = exp(polyval(p, log(x1))) * offset;
    y2 = y1 * (x2 / x1)^p(1);

    args = {'Color', color, 'LineWidth', 1.0, 'HandleVisibility', 'off'};
    plot(ax, [x1 x2], [y1 y2], args{:});    % hypotenuse: the reference slope
    plot(ax, [x1 x2], [y1 y1], args{:});    % horizontal leg, on top
    plot(ax, [x2 x2], [y1 y2], args{:});    % vertical leg, on the right

    text(ax, exp(mean(log([x1 x2]))), y1 * 1.03, sprintf('$%.2f$', p(1)), ...
        'Color', color, 'FontSize', 10, 'HandleVisibility', 'off', ...
        'HorizontalAlignment', 'center', 'VerticalAlignment', 'bottom');

    y_hi = y1 * 1.25;   % include the label in the reported extent
end
