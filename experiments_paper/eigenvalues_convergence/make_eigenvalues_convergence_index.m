function make_eigenvalues_convergence_index(index, name)
%MAKE_EIGENVALUES_CONVERGENCE_INDEX Convergence of one eigenvalue against DOF (paper).
%
%   make_eigenvalues_convergence_index(index) follows the eigenvalue of that
%   index -- lambda_100, say -- of the Dirichlet Laplacian across a sweep of
%   resolutions, computed with DST, finite differences and finite elements, and
%   plots it against the degrees of freedom actually used by each discretisation.
%
%   make_eigenvalues_convergence_index(index, name) restricts the run to the
%   single domain "name". Pass "" or [] to keep all domains.
%
%   Companion of MAKE_EIGENVALUES_CONVERGENCE, which follows the ground state:
%   same domains, same resolution sweeps, same timing and caching. Any index may
%   be asked for, the ground state included, and what changes with the index is
%   the reference. The ground state has a published one, the MPS value of Betcke
%   & Trefethen; deeper in the spectrum nobody has published anything, and the
%   finest DST run of the sweep stands in. REFERENCE_FOR picks between them and
%   the figures say which they used.
%
%   Deep in the spectrum a coarse discretisation has little left to work with:
%   the hundredth eigenvalue of a run carrying a few hundred degrees of freedom
%   sits near the top of its own spectrum, where every method is at its worst.
%   That is the interesting end of the plot, so the sweep is not trimmed at the
%   coarse side; a run whose spectrum is shorter than the index asked for is
%   dropped, and its curve starts at the first resolution that reaches it.
%
%   Output, into results_paper/eigenvalues_convergence/, the index carried in
%   the file names so that several indices can live side by side:
%     - <domain>_convergence_lambda<index>.csv   one row per (method,
%       resolution) with the resolution parameter, the dof count, the eigenvalue
%       and the timing,
%     - <domain>_convergence_lambda<index>.eps   the eigenvalue against DOF on
%       linear axes, one curve per method, with the published value drawn as a
%       level where there is one. No title: the domain is carried by the file
%       name, as elsewhere in experiments_paper, and
%     - <domain>_convergence_lambda<index>_error.eps   the same runs as a
%       relative error against DOF on log-log axes, in the manner of
%       MAKE_EIGENVALUES_CONVERGENCE and with the same slope triangles.
%   Vector EPS, not raster: see SAVE_FIGURE.
%   Against a DST reference the error figure stops the DST curve three runs short
%   of it and gives it no triangle, for the reasons in PLOT_ERROR_RUNS and CURVE.
%   The value figure always shows every run of every method.
%
%   The CSV is a cache, read run by run: a run already in it is taken over as it
%   stands, timing included, and only what is missing is computed. So a figure
%   can be restyled for free, and a resolution added to the sweep costs its own
%   run and nothing else. Delete the CSV to force the whole sweep again.
%
%   Requires the PDE Toolbox (FEM).
%
%   See also MAKE_EIGENVALUES_CONVERGENCE.

    if nargin < 1 || isempty(index)
        error('eigenvalues_convergence:index', ...
            'Give the eigenvalue index, e.g. make_eigenvalues_convergence_index(100).');
    end
    if ~isscalar(index) || index < 1 || index ~= fix(index)
        error('eigenvalues_convergence:index', ...
            'The eigenvalue index must be a positive whole number, got %g.', index);
    end

    here         = fileparts(mfilename('fullpath'));
    project_root = fileparts(fileparts(here));
    run(fullfile(project_root, 'startup.m'));
    addpath(fullfile(project_root, 'experiments'));   % domain_catalog_dst / _fem

    out_dir = fullfile(project_root, 'results_paper', 'eigenvalues_convergence');
    if ~exist(out_dir, 'dir')
        mkdir(out_dir);
    end

    cfgs = domain_configs();
    if nargin >= 2 && ~isempty(name)
        cfgs = cfgs(strcmp({cfgs.name}, name));
        if isempty(cfgs)
            error('eigenvalues_convergence:unknownDomain', ...
                'Unknown domain "%s" (known: L_shaped).', name);
        end
    end
    for i = 1:numel(cfgs)
        cfg = cfgs(i);
        fprintf('=== %s, eigenvalue %d ===\n', cfg.name, index);
        stem = fullfile(out_dir, sprintf('%s_convergence_lambda%d', cfg.name, index));
        csv = [stem '.csv'];
        if exist(csv, 'file')
            cached = read_runs_csv(csv);
            fprintf('  (reading %s: %d runs)\n', csv, numel(cached));
        else
            cached = empty_runs();
        end
        [runs, computed, reused] = compute_runs(cfg, index, cached);
        if computed > 0
            write_runs_csv(csv, cfg, runs, index);
            fprintf('  Wrote %s (%d run(s) computed, %d reused)\n', ...
                    csv, computed, reused);
        end
        ref = reference_for(cfg, index, runs);
        plot_runs(stem, runs, index, ref);
        plot_error_runs([stem '_error'], runs, index, ref);
    end
end


function cfgs = domain_configs()
%DOMAIN_CONFIGS Geometry and resolution sweeps per domain.
%
%   The sweeps follow MAKE_EIGENVALUES_CONVERGENCE, so that the figures of a
%   domain are read against the same DOF counts, with one resolution added at the
%   fine end.
%
%   REF_VALUE is the published eigenvalue, and REF_INDEX the index it belongs to
%   -- the ground state, the only one for which anybody has published a value.
%   Asked for that index, the figures use it: the exact level is drawn, and the
%   error measured against it is a true error, relative to the eigenvalue. Asked
%   for any other, they fall back on the finest DST run of the sweep, with the
%   caveats in PLOT_ERROR_RUNS.
    cfgs = struct('name', {}, 'pretty', {}, 'box', {}, 'phi', {}, ...
                  'gd', {}, 'ns', {}, 'sf', {}, 'M_grid', {}, 'Hmax_fem', {}, ...
                  'ref_index', {}, 'ref_value', {}, 'ref_label', {});

    % --- L-shaped domain ----------------------------------------------------
    % Box [-1,1]^2, so h = 2/(M+1) and M+1 must be even for the re-entrant
    % corner at the origin to fall on a grid line -- alignment matters a lot
    % here, the corner carrying the leading singularity. The M values sweep the
    % DOF counts from about 100 to about 15000: steps of roughly 500 up to 5000,
    % then coarser ones, the dense eig cost growing as dofs^3. The Hmax values
    % target the same counts through the empirical dofs ~ 12.8/Hmax^2.
    %
    % The sweep runs one resolution past the ground-state one, to M = 141. The
    % finest DST run stands in for the reference value that does not exist away
    % from the ground state, so it has to sit a good way ahead of the runs
    % measured against it. On this grid dofs = 3*M^2/4 - M/2 exactly, so M = 141
    % carries 14840 of them against 9861 at M = 115, and costs about fourteen
    % minutes on its own -- roughly as much as the whole rest of the sweep.
    %
    % The FEM meshes follow the same DOF counts, through the empirical
    % dofs*Hmax^2 ~ 13.6 that the finer meshes here obey to within a few per
    % cent, and the last one matches the finest DST run: Hmax = 0.030 for about
    % 15000 degrees of freedom.
    cfgs(end+1) = struct( ...
        'name', 'L_shaped', 'pretty', 'L-shaped', ...
        'box', [-1 1 -1 1], ...
        'phi', @(x, y) indicator_L_shaped(x, y, -1, 0, 1, -1, 0, 1), ...
        'gd', [[3; 4; -1; 1; 1; -1; -1; -1; 1; 1], [3; 4; 0; 1; 1; 0; 0; 0; 1; 1]], ...
        'ns', char('R1', 'R2')', 'sf', 'R1-R2', ...
        'M_grid',   [13 17 31 41 49 55 63 69 75 81 93 105 115 141], ...
        'Hmax_fem', [0.33 0.25 0.135 0.10 0.085 0.0755 0.066 0.060 0.055 0.051 ...
                     0.0446 0.0395 0.036 0.030], ...
        'ref_index', 1, 'ref_value', 9.6397238440219, 'ref_label', 'MPS');
end


function [runs, computed, reused] = compute_runs(cfg, index, cached)
%COMPUTE_RUNS One timed eigenvalue per method and resolution, cached run by run.
%
%   CACHED holds whatever the CSV already had. A planned run found there is
%   taken over as it stands, timing included; only what is missing is computed.
%   That is what makes the sweep extensible: a resolution added to the config
%   costs its own run and nothing else, where an all-or-nothing cache would put
%   the whole half hour back on the bill. Runs are returned in the order of the
%   plan, so that the rewritten CSV keeps the shape of the old one.
%
%   COMPUTED counts the computations performed and REUSED the runs taken from the
%   cache. They need not add up to the table: a run too coarse to reach the index
%   is computed and then dropped. COMPUTED is zero when the CSV already covered
%   the whole sweep, and the caller then leaves the file alone.
%
%   A run too coarse to reach the index never enters the CSV, so it is attempted
%   again on each pass. Those are the cheapest runs of the sweep, a fraction of a
%   second, and the alternative is a cache that records absences.
%
%   Warm-up (see WARM_UP) is paid only if something is to be computed, so that a
%   figure-only regeneration stays free.
    [x_range, y_range] = bounding_box(cfg.box(1), cfg.box(2), cfg.box(3), cfg.box(4));

    plan = build_plan(cfg);
    hits = arrayfun(@(p) find_run(cached, p.method, p.resolution), plan);
    if all(hits)
        fprintf('  every run of the sweep is in the CSV\n');
    else
        warm_up(cfg, x_range, y_range);
    end

    runs = empty_runs();
    computed = 0;
    reused = 0;
    for i = 1:numel(plan)
        p = plan(i);
        if hits(i)
            runs(end+1) = pick_run(cached, p.method, p.resolution); %#ok<AGROW>
            reused = reused + 1;
            continue;
        end
        t = tic;
        [evals, dofs] = p.compute(x_range, y_range);
        tsec = toc(t);
        computed = computed + 1;
        if numel(evals) < index
            fprintf('  %-4s %-12s dofs = %-5d  only %d eigenvalues, skipped\n', ...
                    p.method, p.resolution, dofs, numel(evals));
            continue;
        end
        runs(end+1) = mk(p.method, p.resolution, dofs, evals(index), tsec); %#ok<AGROW>
        report(runs(end), index);
    end
end


function plan = build_plan(cfg)
%BUILD_PLAN Every run of the sweep, in output order, with how to compute it.
%
%   DST and FD walk M_grid together, FEM walks Hmax_fem after them. The
%   resolution string is the cache key, so it has to be written the same way here
%   and in the CSV -- hence one place that builds it.
    plan = struct('method', {}, 'resolution', {}, 'compute', {});
    for k = 1:numel(cfg.M_grid)
        M = cfg.M_grid(k);
        res = sprintf('M = %d', M);
        plan(end+1) = struct('method', 'DST', 'resolution', res, ...
            'compute', @(xr, yr) dst_spectrum(xr, yr, M, cfg.phi)); %#ok<AGROW>
        plan(end+1) = struct('method', 'FD', 'resolution', res, ...
            'compute', @(xr, yr) fd_spectrum(cfg, M)); %#ok<AGROW>
    end
    for k = 1:numel(cfg.Hmax_fem)
        h = cfg.Hmax_fem(k);
        plan(end+1) = struct('method', 'FEM', 'resolution', sprintf('Hmax = %g', h), ...
            'compute', @(xr, yr) fem_spectrum(cfg, h)); %#ok<AGROW>
    end
end


function [evals, dofs] = dst_spectrum(x_range, y_range, M, phi)
    [L, info] = make_dst_laplace_mat_batched(x_range, y_range, M, phi);
    evals = sort(-real(eig(L)), 'ascend');
    dofs = info.dofs;
end


function [evals, dofs] = fd_spectrum(cfg, M)
    [evals, info] = fd_laplace_spectrum(struct('box', cfg.box, 'phi', cfg.phi, 'M', M));
    evals = sort(real(evals(:)), 'ascend');
    dofs = info.dofs;
end


function [evals, dofs] = fem_spectrum(cfg, h)
    entry = struct('gd', cfg.gd, 'ns', cfg.ns, 'sf', cfg.sf, ...
                   'Hmax_eig', h, 'Hmax_solvepdeeig', h);
    [evals, info] = fem_laplace_spectrum(entry, "eig");
    evals = sort(real(evals(:)), 'ascend');
    dofs = info.dofs;
end


function tf = find_run(runs, method, resolution)
%FIND_RUN Whether the cache holds this run.
    tf = ~isempty(runs) && any(strcmp({runs.method}, method) & ...
                               strcmp({runs.resolution}, resolution));
end


function r = pick_run(runs, method, resolution)
%PICK_RUN The cached run, the first if a hand-edited CSV repeats one.
    hit = find(strcmp({runs.method}, method) & strcmp({runs.resolution}, resolution), 1);
    r = runs(hit);
end


function runs = empty_runs()
    runs = struct('method', {}, 'resolution', {}, 'dofs', {}, 'lambda', {}, 'time', {});
end


function warm_up(cfg, x_range, y_range)
%WARM_UP Discarded runs of each method on this domain, at its cheapest resolution.
%
%   Warm-up has two components, both measured on this code: MATLAB's toolbox
%   load and JIT compilation, paid once per session (about 1.7 s for FEM), and a
%   smaller per-geometry cost, paid again by the first call on each new domain.
%   Left in the sweep, they land on the coarsest run and put the smallest DOF
%   count at the top of the timings. Two calls rather than one, because the
%   first geometry of a session is still decaying on its second call.
%
%   Eager rather than lazy here: COMPUTE_RUNS is reached only when the CSV is
%   missing, that is when every run is recomputed anyway.
    fprintf('  warm-up: DST, FD, FEM at the coarsest resolution (discarded)\n');
    M = cfg.M_grid(1);
    h = cfg.Hmax_fem(1);
    entry = struct('gd', cfg.gd, 'ns', cfg.ns, 'sf', cfg.sf, ...
                   'Hmax_eig', h, 'Hmax_solvepdeeig', h);
    for i = 1:2
        L = make_dst_laplace_mat_batched(x_range, y_range, M, cfg.phi);
        eig(L);
        fd_laplace_spectrum(struct('box', cfg.box, 'phi', cfg.phi, 'M', M));
        fem_laplace_spectrum(entry, "eig");
    end
end


function s = mk(method, resolution, dofs, lambda, tsec)
    s = struct('method', method, 'resolution', resolution, 'dofs', dofs, ...
               'lambda', lambda, 'time', tsec);
end


function report(r, index)
    fprintf('  %-4s %-12s dofs = %-5d  lambda_%d = %.8f  time = %6.2f s\n', ...
            r.method, r.resolution, r.dofs, index, r.lambda, r.time);
end


function write_runs_csv(csv, cfg, runs, index)
%WRITE_RUNS_CSV One row per run: method, resolution, dofs, eigenvalue, time.
    fid = fopen(csv, 'w');
    if fid == -1
        error('eigenvalues_convergence:csv', 'Could not open %s for writing.', csv);
    end
    closer = onCleanup(@() fclose(fid));
    fprintf(fid, ['# Domain: %s (convergence of eigenvalue %d against DOF, ', ...
                  'dense eig)\n'], cfg.name, index);
    fprintf(fid, '# No reference value: the figure plots lambda_%d itself.\n', index);
    fprintf(fid, 'method,resolution,dofs,lambda_%d,time_s\n', index);
    for i = 1:numel(runs)
        fprintf(fid, '%s,%s,%d,%.12g,%.4f\n', runs(i).method, runs(i).resolution, ...
                runs(i).dofs, runs(i).lambda, runs(i).time);
    end
end


function runs = read_runs_csv(csv)
%READ_RUNS_CSV Read back a run table written by WRITE_RUNS_CSV.
    fid = fopen(csv, 'r');
    if fid == -1
        error('eigenvalues_convergence:csvread', 'Could not open %s.', csv);
    end
    closer = onCleanup(@() fclose(fid));
    runs = empty_runs();
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


function ref = reference_for(cfg, index, runs)
%REFERENCE_FOR The value the figures measure against, and what it is worth.
%
%   PUBLISHED, when the index asked for is the one the config carries a value
%   for: the ground state, against the MPS value of Betcke & Trefethen. It is
%   external to the sweep, so every run can be measured against it and the error
%   is a true error.
%
%   Otherwise the finest DST run of the sweep. That is not a true reference, and
%   what it gives is a difference between two computed numbers rather than a
%   distance from the exact value, so the DST runs nearest it are held back from
%   the error figure.
%
%   Either way the error is reported relative to the reference. The eigenvalues
%   run from 9.6 at the ground state to 919 deep in the spectrum, and an absolute
%   error carries that scale with it; dividing it out is what lets the figures of
%   different indices be read against each other.
    if ~isempty(cfg.ref_index) && index == cfg.ref_index
        ref = struct('value', cfg.ref_value, 'label', cfg.ref_label, ...
                     'published', true, 'relative', true, 'dst_trim', 0);
        fprintf('  reference: %s, lambda_%d = %.13f\n', ref.label, index, ref.value);
        return;
    end
    dst = strcmp({runs.method}, 'DST');
    if ~any(dst)
        error('eigenvalues_convergence:noReference', ...
            'No published value for lambda_%d and no DST run to stand in.', index);
    end
    dst_runs = runs(dst);
    [ref_dofs, imax] = max([dst_runs.dofs]);
    ref = struct('value', dst_runs(imax).lambda, 'label', 'DST', ...
                 'published', false, 'relative', true, 'dst_trim', 3);
    fprintf('  reference: DST at %d dofs, lambda_%d = %.9f\n', ...
            ref_dofs, index, ref.value);
end


function plot_runs(stem, runs, index, ref)
%PLOT_RUNS The eigenvalue against DOF, one curve per method.
%
%   The computed eigenvalue itself on linear axes, not an error and not a log-log
%   rate plot: the question the figure answers is whether the three methods
%   settle on a common level and how soon. Where that level is known -- the
%   ground state, see REFERENCE_FOR -- it is drawn as a horizontal line, and the
%   curves are read against it rather than against each other. Colours and
%   markers are those of MAKE_EIGENVALUES_CONVERGENCE, so that the figures of a
%   domain read as one family.
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

    for mi = 1:numel(methods)
        [d, v] = curve(runs, methods{mi});
        plot(ax, d, v, ['-' markers{mi}], 'Color', colors{mi}, 'LineWidth', 2.0, ...
            'MarkerSize', 7, 'MarkerFaceColor', 'w', ...
            'DisplayName', sprintf('\\texttt{%s}', methods{mi}));
    end
    if ref.published
        % Drawn last so that it lies over the curves, and dashed and black so
        % that it reads as the level they are converging to and not as a fourth
        % method.
        yline(ax, ref.value, 'k--', 'LineWidth', 1.5, ...
            'DisplayName', sprintf('\\texttt{%s}', ref.label));
    end

    hold(ax, 'off');
    xlabel(ax, 'degrees of freedom');
    ylabel(ax, sprintf('$\\lambda_{%d}$', index));
    % No title: the domain is identified by the output file name.
    legend(ax, 'Location', 'northeast', 'FontSize', 10, 'Box', 'off');
    grid(ax, 'on'); box(ax, 'on');
    set(ax, 'FontSize', 12);

    save_figure(fig, stem);
end


function plot_error_runs(stem, runs, index, ref)
%PLOT_ERROR_RUNS Error against DOF on log-log axes, in the manner of the sibling.
%
%   The relative error, always, so that the figures of different indices can be
%   read against each other; REFERENCE_FOR says what it is relative to and how
%   much that reference is worth.
%
%   Two cases, told apart there. With the published value -- the ground state,
%   the MPS value of Betcke & Trefethen -- the reference lies outside the sweep,
%   so every run can be measured against it over the whole length of its curve
%   and each method carries a slope triangle.
%
%   Without one the finest DST run stands in. The DST curve then stops three runs
%   short of it -- see CURVE -- a point measured against a neighbouring
%   resolution saying how fast DST is still moving rather than how far it is from
%   the true eigenvalue, and it carries no triangle: away from the ground state
%   DST converges before the sweep starts, and a line fitted to a curve that
%   drops once and then lies flat describes the drop, not the method. FD and FEM
%   are drawn whole either way.
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

    % Each triangle sits above its own curve, staggered so that no two stack up.
    % With the published reference all three curves are power laws and all three
    % get one, in the placement MAKE_EIGENVALUES_CONVERGENCE settled on for this
    % data: FD takes the coarse left end where its gap to DST is still wide, DST
    % the middle, FEM the fine right end where it has pulled away from FD.
    % Against a DST reference, DST gets none.
    if ref.published
        triangle = [true, true, true];
        % FEM keeps a shorter span than the other two. It is the lowest curve and
        % the steepest, so a triangle drawn over the usual quarter of the range
        % grows taller than the gap to FD above it and its top edge cuts across
        % that curve; the height goes with the span, so the span comes in.
        spans    = {[0.32 0.55], [0.08 0.28], [0.74 0.90]};
        offsets  = [1.13, 1.25, 1.10];
    else
        triangle = [false, true, true];
        spans    = {[], [0.38 0.61], [0.68 0.91]};
        offsets  = [NaN, 2.2, 1.8];
    end

    lo = inf; hi = 0;
    for mi = 1:numel(methods)
        [d, v] = curve(runs, methods{mi}, ref.dst_trim);
        err = abs(v - ref.value);
        if ref.relative
            err = err / abs(ref.value);
        end
        keep = err > 0;          % belt and braces: no zero can reach a log axis
        d = d(keep); err = err(keep);
        loglog(ax, d, err, ['-' markers{mi}], 'Color', colors{mi}, 'LineWidth', 2.0, ...
            'MarkerSize', 7, 'MarkerFaceColor', 'w', ...
            'DisplayName', sprintf('\\texttt{%s}', methods{mi}));
        yhi = 0;
        if triangle(mi)
            % Least-squares algebraic rate: err ~ dofs^p, annotated by the triangle.
            p = polyfit(log(d), log(err), 1);
            yhi = slope_triangle(ax, d, p, colors{mi}, spans{mi}, offsets(mi));
        end
        lo = min([lo, err]);
        hi = max([hi, err, yhi]);
    end

    hold(ax, 'off');
    set(ax, 'XScale', 'log', 'YScale', 'log');
    ylim(ax, [lo * 0.72, hi * 1.35]);
    xlabel(ax, 'degrees of freedom');
    if ref.relative
        % A built fraction rather than a solidus, as in the sub-captions of the
        % snippets. It is set a size larger than the other axis text, the
        % numerator and denominator of a fraction being rendered small.
        ylabel(ax, sprintf(['$\\frac{|\\lambda_{%d} - \\lambda_{%d}^{\\mathrm{%s}}|}', ...
                            '{\\lambda_{%d}^{\\mathrm{%s}}}$'], ...
                           index, index, ref.label, index, ref.label), 'FontSize', 16);
    else
        ylabel(ax, sprintf('$|\\lambda_{%d} - \\lambda_{%d}^{\\mathrm{%s}}|$', ...
                           index, index, ref.label));
    end
    % No title: the domain is identified by the output file name.
    legend(ax, 'Location', 'southwest', 'FontSize', 10, 'Box', 'off');
    grid(ax, 'on'); box(ax, 'on');
    set(ax, 'FontSize', 12);

    save_figure(fig, stem);
end


function save_figure(fig, stem)
%SAVE_FIGURE One figure as vector EPS, then close it.
%
%   These plots are line art -- axes, curves, markers, a few labels -- so a
%   vector figure is both smaller than a raster of the same page area and sharp
%   at any size, and there is no resolution to choose. Read by latex/dvips
%   directly, and by pdflatex through epstopdf.
    eps_file = [stem '.eps'];
    try
        exportgraphics(fig, eps_file, 'ContentType', 'vector');
    catch
        print(fig, eps_file, '-depsc2', '-painters');
    end
    fprintf('  Wrote %s\n', eps_file);
    close(fig);
end


function [d, v] = curve(runs, method, dst_trim)
%CURVE Sorted DOF and eigenvalue of one method, ready to plot.
%
%   DST_TRIM, none by default, holds back that many of the finest DST runs. Only
%   the error figure asks for it, its reference being a DST run itself: the
%   finest would sit at zero error, and the two below it are close enough to it
%   in resolution that their difference says how fast DST is still moving rather
%   than how far it is from the true eigenvalue. The value figure takes the whole
%   sweep, every run there being just a computed eigenvalue, with no reference
%   for it to be measured against. Either way the CSV keeps everything -- what is
%   trimmed is the drawing, not the data.
    if nargin < 3
        dst_trim = 0;
    end

    sel = strcmp({runs.method}, method);
    d = [runs(sel).dofs];
    v = [runs(sel).lambda];
    [d, ord] = sort(d);
    v = v(ord);
    if strcmp(method, 'DST') && dst_trim > 0
        n = max(0, numel(d) - dst_trim);
        d = d(1:n);
        v = v(1:n);
    end
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
