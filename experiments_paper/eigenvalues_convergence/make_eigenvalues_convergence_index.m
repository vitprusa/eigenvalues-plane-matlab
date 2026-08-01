function make_eigenvalues_convergence_index(index, name)
%MAKE_EIGENVALUES_CONVERGENCE_INDEX Convergence of one eigenvalue against DOF (paper).
%
%   make_eigenvalues_convergence_index(index) follows the eigenvalue of that
%   index -- lambda_100, say -- of the Dirichlet Laplacian across a sweep of
%   resolutions, computed with DST, finite differences and finite elements, and
%   plots it against the degrees of freedom actually used by each discretisation.
%
%   make_eigenvalues_convergence_index(index, name) restricts the run to the
%   single domain "name" -- rectangle, isosceles_triangle or L_shaped. Pass ""
%   or [] to keep all three.
%
%   Any index may be asked for, the ground state included, and what changes with
%   the domain and the index is the reference. The rectangle and the right
%   isosceles triangle have a closed-form spectrum, so every index has an exact
%   reference. The L-shape has a published value for the ground state alone, the
%   MPS value of Betcke & Trefethen; deeper in its spectrum nobody has published
%   anything, and the finest DST run of the sweep stands in. REFERENCE_FOR picks
%   between the three and the figures say which they used.
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
%       MAKE_EIGENVALUES_CONVERGENCE and with the same slope triangles, and
%     - spectra/<domain>_<method>_<resolution>-eigenvalues.csv   the whole
%       spectrum of each run, which is the cache the other three are built from.
%   Vector EPS, not raster: see SAVE_FIGURE.
%   Against a DST reference the error figure gives the DST curve no triangle, for
%   the reason in PLOT_ERROR_RUNS. The whole of every curve is drawn either way:
%   the reference is a run of its own, finer than the sweep, not a member of it.
%
%   On the rectangle the DST error is at roundoff -- the domain fills its
%   bounding box, so the sine basis the operator is built in is the exact
%   eigenbasis -- and that curve gets no triangle either: see PLOT_ERROR_RUNS.
%
%   The cache is the whole spectrum of each run, one file per run under
%   results_paper/eigenvalues_convergence/spectra/ (see READ_SPECTRUM). What
%   costs the time here is the dense eig, and it returns the whole spectrum
%   whatever index is asked for; keeping only one eigenvalue of it made every
%   further index pay for the sweep again. Kept, the sweep is paid once per
%   domain and every index after the first is a redraw, at any index whatever.
%   Delete the spectra of a domain to force it again.
%
%   A run whose spectrum is shorter than the index asked for is left out of the
%   figures, and its stored spectrum says so, so it is not attempted again on the
%   next pass.
%
%   The per-index run tables written next to the figures are output, not cache.
%   The exception is the legacy path: where no spectrum has been stored but an
%   older run table holds the run, its eigenvalue and timing are taken from
%   there, which is what keeps the L-shape's sweep -- half an hour of it, already
%   computed and committed -- from having to be run again.
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
    % The spectra are the cache, and there is one file per run of every sweep:
    % their own folder, so that the figures and the run tables stay legible next
    % to each other in the parent.
    spectra_dir = fullfile(out_dir, 'spectra');
    if ~exist(spectra_dir, 'dir')
        mkdir(spectra_dir);
    end

    cfgs = domain_configs();
    if nargin >= 2 && ~isempty(name)
        known = {cfgs.name};
        cfgs = cfgs(strcmp(known, name));
        if isempty(cfgs)
            error('eigenvalues_convergence:unknownDomain', ...
                'Unknown domain "%s" (known: %s).', name, strjoin(known, ', '));
        end
    end
    for i = 1:numel(cfgs)
        cfg = cfgs(i);
        fprintf('=== %s, eigenvalue %d ===\n', cfg.name, index);
        stem = fullfile(out_dir, sprintf('%s_convergence_lambda%d', cfg.name, index));
        csv = [stem '.csv'];
        if exist(csv, 'file')
            legacy = read_runs_csv(csv);
        else
            legacy = empty_runs();
        end
        [runs, computed, reused] = compute_runs(cfg, index, spectra_dir, legacy);
        ref = reference_for(cfg, index, dst_reference(cfg, index, spectra_dir));
        write_runs_csv(csv, cfg, runs, index, ref);
        fprintf('  Wrote %s (%d run(s) computed, %d from cache)\n', ...
                csv, computed, reused);
        plot_runs(stem, runs, index, ref);
        plot_error_runs([stem '_error'], runs, index, ref);
    end
end


function cfgs = domain_configs()
%DOMAIN_CONFIGS Geometry and resolution sweeps per domain.
%
%   Every sweep covers the same band of DOF counts, from about 100 to about
%   10000, so that the figures of the three domains are read against each other;
%   the L-shape and the four catalogue domains run one resolution past it, to
%   about 15000. Their reference is a separate run again, at M_REF; see
%   DST_REFERENCE.
%
%   EXACT_FUN, where the spectrum is known in closed form, returns the exact
%   eigenvalue of the index asked for. It is external to the sweep and available
%   at every index, so on those domains the figures always draw the exact level
%   and always measure a true error against it.
%
%   Without one, REF_VALUE is the published eigenvalue and REF_INDEX the index it
%   belongs to -- the ground state, the only one for which anybody has published
%   a value on the L-shape. Asked for that index, the figures use it. Asked for
%   any other, they fall back on the finest DST run of the sweep, with the
%   caveats in PLOT_ERROR_RUNS.
    cfgs = struct('name', {}, 'pretty', {}, 'box', {}, 'phi', {}, ...
                  'gd', {}, 'ns', {}, 'sf', {}, 'M_grid', {}, 'Hmax_fem', {}, ...
                  'M_ref', {}, 'exact_fun', {}, 'ref_index', {}, 'ref_value', {}, ...
                  'ref_label', {});

    % --- rectangle [0, 2*pi] x [0, pi] ---------------------------------------
    % Box [0,2pi]x[0,pi], so h = 2*pi/(M+1) and M+1 must be even for the top edge
    % y = pi to fall on a grid line. The domain then fills its bounding box and
    % carries M columns of (M-1)/2 interior points, dofs = M*(M-1)/2 exactly, so
    % the M values below sweep 105 to 10731 degrees of freedom.
    %
    % There is no run ahead of the sweep here: the reference is the closed-form
    % eigenvalue, m^2/4 + n^2, so nothing has to stand in for it and the sweep
    % stops where the L-shape's measured curves stop.
    %
    % Filling the bounding box is what makes this domain the control of the set.
    % The DST operator is built in the sine basis of the box, which is the exact
    % eigenbasis of the rectangle, so DST returns m^2/4 + n^2 to roundoff at every
    % resolution -- what it can get wrong is only which modes the grid carries,
    % and that shows up at indices deep enough to reach the edge of the grid.
    %
    % The FEM meshes follow the same DOF counts through the empirical
    % dofs*Hmax^2 ~ 88 that the finer meshes here obey to within a few per cent.
    cfgs(end+1) = struct( ...
        'name', 'rectangle', 'pretty', 'rectangular', ...
        'box', [0 2*pi 0 pi], ...
        'phi', @(x, y) indicator_rectangle(x, y, 0, 2*pi, 0, pi), ...
        'gd', [3; 4; 0; 2*pi; 2*pi; 0; 0; 0; pi; pi], ...
        'ns', char('R1')', 'sf', 'R1', ...
        'M_grid',   [15 21 31 41 51 61 71 81 91 105 119 133 147], ...
        'Hmax_fem', [0.85 0.61 0.42 0.32 0.26 0.22 0.188 0.165 0.147 0.127 ...
                     0.1125 0.100 0.0914], ...
        'M_ref', [], 'exact_fun', @exact_rectangle, ...
        'ref_index', [], 'ref_value', [], 'ref_label', '');

    % --- right isosceles triangle, legs pi -----------------------------------
    % Box [0,pi]^2, so h = pi/(M+1), and the hypotenuse y = x runs through grid
    % points whatever M is: no alignment condition to meet. The interior points
    % are those with y < x, dofs = M*(M-1)/2 exactly, the same counts as the
    % rectangle, and the same M values sweep them.
    %
    % Reference: lambda = m^2 + n^2 with m > n >= 1, again exact at every index,
    % so again no run ahead of the sweep.
    %
    % The FEM meshes follow the DOF counts through the empirical dofs*Hmax^2 ~ 22
    % of this geometry -- a quarter of the rectangle's, the domain being a quarter
    % of its area -- obeyed to within a few per cent by the finer meshes.
    cfgs(end+1) = struct( ...
        'name', 'isosceles_triangle', 'pretty', 'right isosceles triangle', ...
        'box', [0 pi 0 pi], ...
        'phi', @(x, y) indicator_isosceles_triangle(x, y, 0, pi, 0), ...
        'gd', [2; 3; 0; pi; pi; 0; 0; pi], ...
        'ns', char('T1')', 'sf', 'T1', ...
        'M_grid',   [15 21 31 41 51 61 71 81 91 105 119 133 147], ...
        'Hmax_fem', [0.43 0.30 0.21 0.16 0.13 0.109 0.094 0.082 0.073 0.063 ...
                     0.056 0.050 0.0456], ...
        'M_ref', [], 'exact_fun', @exact_isosceles_triangle, ...
        'ref_index', [], 'ref_value', [], 'ref_label', '');

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
        'M_ref', 163, 'exact_fun', [], ...
        'ref_index', 1, 'ref_value', 9.6397238440219, 'ref_label', 'MPS');

    % --- the domains with no closed form and no analytic geometry to spell out -
    % Ellipse minus a quadrant, H, and the two GWW isospectral drums. Their box
    % and indicator come from DOMAIN_CATALOG_DST and their decsg geometry from
    % DOMAIN_CATALOG_FEM, as in MAKE_EIGENVALUES_HEAD_PAPER_TABLES, rather than
    % being written out again here.
    %
    % Like the L-shape they have no closed-form spectrum, so the ground state is
    % measured against the published MPS value of Betcke & Trefethen kept in
    % data/, and every other index against the finest DST run of the sweep. That
    % is why each sweep carries a fourteenth resolution past the band the other
    % domains stop at: the run that stands in for a reference has to sit well
    % ahead of the runs measured against it.
    %
    % The M values keep the domain edges on grid lines, which is a different
    % condition on each: the grid spacing is the width of the bounding box over
    % M+1, and an edge at distance d from its left needs d/h whole. They sweep
    % the same band of DOF counts as the other domains, from about 100 to about
    % 10500, and one run beyond.
    %
    % The Hmax values target those same DOF counts through the empirical
    % dofs*Hmax^2 of each geometry -- 21.4 for the ellipse, 31 for H, 62.5 for
    % the drums -- which the finer meshes obey to within a few per cent. The
    % coarsest two or three of each run a little under it and are set from a
    % measured mesh rather than from the constant.
    %
    % Columns: catalog name for the two lookups, output name (file names and the
    % -name filter; H is H_shaped there, as in results_paper/eigenvalues_head/),
    % caption name, DST/FD resolutions, FEM mesh sizes, published lambda_1, and
    % the resolution of the reference run (see DST_REFERENCE).
    extra = { ...
        'ellipse_minus_quadrant', 'ellipse_minus_quadrant', 'ellipse-minus-quadrant', ...
            [19 27 43 59 71 83 95 107 119 135 151 171 191 223], ...
            [0.40 0.30 0.190 0.143 0.119 0.1025 0.0897 0.0797 0.072 0.0635 ...
             0.057 0.0504 0.0452 0.0387], 5.868746216295, 259; ...
        'H', 'H_shaped', 'H-shaped', ...
            [14 17 29 38 47 53 62 68 74 83 92 101 116 137], ...
            [0.40 0.36 0.21 0.1675 0.128 0.12 0.103 0.094 0.086 0.077 ...
             0.0678 0.0632 0.055 0.0466], 7.7330888559, 161; ...
        'gww1', 'gww1', 'GWW1 isospectral drum', ...
            [17 23 35 47 59 65 77 89 101 113 125 137 155 179], ...
            [0.71 0.54 0.355 0.262 0.21 0.194 0.1667 0.144 0.1267 0.113 ...
             0.102 0.093 0.082 0.0712], 2.537943999798, 227; ...
        'gww2', 'gww2', 'GWW2 isospectral drum', ...
            [17 23 35 47 59 65 77 89 101 113 125 137 155 179], ...
            [0.71 0.54 0.355 0.262 0.21 0.194 0.1667 0.144 0.1267 0.113 ...
             0.102 0.093 0.082 0.0712], 2.537943999798, 227 ...
    };
    for i = 1:size(extra, 1)
        dstc = domain_catalog_dst(extra{i, 1});
        femc = domain_catalog_fem(extra{i, 1});
        cfgs(end+1) = struct( ...
            'name', extra{i, 2}, 'pretty', extra{i, 3}, ...
            'box', dstc.box, 'phi', dstc.phi, ...
            'gd', femc.gd, 'ns', femc.ns, 'sf', femc.sf, ...
            'M_grid', extra{i, 4}, 'Hmax_fem', extra{i, 5}, ...
            'M_ref', extra{i, 7}, 'exact_fun', [], ...
            'ref_index', 1, 'ref_value', extra{i, 6}, 'ref_label', 'MPS'); %#ok<AGROW>
    end
end


function lambda = exact_rectangle(index)
%EXACT_RECTANGLE The index-th eigenvalue of [0,2pi]x[0,pi], in closed form.
%
%   lambda_{m,n} = m^2/4 + n^2, m, n >= 1, sorted ascending with multiplicity.
%   The pairs are enumerated below a bound rather than over a fixed rectangle of
%   (m, n): a truncation in m and n alone can miss an eigenvalue smaller than one
%   it keeps, and then the index-th entry of the sorted list is not the index-th
%   eigenvalue. Below a bound nothing is missed, and the bound is doubled until
%   it holds enough of the spectrum.
    bound = 64;
    while true
        [m, n] = meshgrid(1:ceil(2 * sqrt(bound)), 1:ceil(sqrt(bound)));
        v = m.^2 / 4 + n.^2;
        v = sort(v(v <= bound));
        if numel(v) >= index
            lambda = v(index);
            return;
        end
        bound = 2 * bound;
    end
end


function lambda = exact_isosceles_triangle(index)
%EXACT_ISOSCELES_TRIANGLE The index-th eigenvalue of the triangle, in closed form.
%
%   Right isosceles triangle with legs pi: lambda_{m,n} = m^2 + n^2 with
%   m > n >= 1, the antisymmetric half of the square's spectrum. Enumerated below
%   a doubling bound, for the reason in EXACT_RECTANGLE.
    bound = 64;
    while true
        [m, n] = meshgrid(1:ceil(sqrt(bound)), 1:ceil(sqrt(bound)));
        v = m.^2 + n.^2;
        v = sort(v(m > n & v <= bound));
        if numel(v) >= index
            lambda = v(index);
            return;
        end
        bound = 2 * bound;
    end
end


function [runs, computed, reused] = compute_runs(cfg, index, spectra_dir, legacy)
%COMPUTE_RUNS One timed eigenvalue per method and resolution, off the spectra.
%
%   Three places a run can come from, in this order:
%
%   The stored spectrum, which is complete, so it either holds the index or the
%   run has no such eigenvalue. The whole point of the restructuring: the dense
%   eig returns the spectrum, so the spectrum is what is kept, and the second
%   index asked for on a domain costs nothing.
%
%   LEGACY, the run table of this index from before there were spectra, for the
%   runs the spectra do not hold. It carries one eigenvalue and its timing, which
%   is all this index needs, and it saves recomputing the L-shape.
%
%   The computation, otherwise. What comes back is stored as a spectrum, so it is
%   paid once.
%
%   COMPUTED counts the computations performed and REUSED the runs taken from
%   either cache. They need not add up to the table: a run whose spectrum is
%   shorter than the index is left out of it.
%
%   Warm-up (see WARM_UP) is paid only if something is to be computed, so that a
%   figure-only regeneration stays free.
%
%   Runs are returned in the order of the plan, so that the rewritten run table
%   keeps the shape of the old one.
    [x_range, y_range] = bounding_box(cfg.box(1), cfg.box(2), cfg.box(3), cfg.box(4));

    plan = build_plan(cfg);
    paths = arrayfun(@(p) spectrum_path(spectra_dir, cfg, p.method, p.resolution), ...
                     plan, 'UniformOutput', false);
    from_cache = arrayfun(@(i) ...
        spectrum_serves(paths{i}, index) || find_run(legacy, plan(i).method, plan(i).resolution), ...
        1:numel(plan));
    if all(from_cache)
        fprintf('  every run of the sweep is cached\n');
    else
        warm_up(cfg, x_range, y_range);
    end

    runs = empty_runs();
    computed = 0;
    reused = 0;
    for i = 1:numel(plan)
        p = plan(i);
        [lambda, dofs, tsec, status] = read_spectrum(paths{i}, index);
        switch status
            case 'hit'
                runs(end+1) = mk(p.method, p.resolution, dofs, lambda, tsec); %#ok<AGROW>
                reused = reused + 1;
                report(runs(end), index, '(cached)');
                continue;
            case 'short'
                fprintf('  %-4s %-12s dofs = %-5d  only %d eigenvalues, left out\n', ...
                        p.method, p.resolution, dofs, lambda);   % lambda: the count
                reused = reused + 1;
                continue;
        end
        if find_run(legacy, p.method, p.resolution)
            runs(end+1) = pick_run(legacy, p.method, p.resolution); %#ok<AGROW>
            reused = reused + 1;
            report(runs(end), index, '(run table)');
            continue;
        end
        t = tic;
        [evals, dofs] = p.compute(x_range, y_range);
        tsec = toc(t);
        computed = computed + 1;
        write_spectrum(paths{i}, cfg, p, dofs, tsec, evals);
        if numel(evals) < index
            fprintf('  %-4s %-12s dofs = %-5d  only %d eigenvalues, left out\n', ...
                    p.method, p.resolution, dofs, numel(evals));
            continue;
        end
        runs(end+1) = mk(p.method, p.resolution, dofs, evals(index), tsec); %#ok<AGROW>
        report(runs(end), index, '');
    end
end


function path = spectrum_path(spectra_dir, cfg, method, resolution)
%SPECTRUM_PATH Where the spectrum of one run of one domain is kept.
%
%   The resolution goes into the name as it is written everywhere else, with the
%   spaces and the equals sign taken out: <domain>_<method>_M147, and
%   <domain>_<method>_Hmax0.0914 for the meshes. One run per file, as in
%   results_paper/eigenvalues_head/.
    slug = regexprep(resolution, '\s*=\s*', '');
    path = fullfile(spectra_dir, sprintf('%s_%s_%s-eigenvalues.csv', ...
                                         cfg.name, lower(method), slug));
end


function tf = spectrum_serves(path, index)
%SPECTRUM_SERVES Whether the stored spectrum settles this index without computing.
    [~, ~, ~, status] = read_spectrum(path, index);
    tf = ~strcmp(status, 'miss');
end


function [lambda, dofs, tsec, status] = read_spectrum(path, index)
%READ_SPECTRUM The index-th eigenvalue of a stored run, if it is there.
%
%   STATUS is "hit" when the spectrum holds the index, and LAMBDA is the
%   eigenvalue; "short" when the run has fewer eigenvalues than the index --
%   LAMBDA is then the number of eigenvalues the run has, for the message -- and
%   "miss" when there is no usable file, which is the one case that has to be
%   computed again.
%
%   A file is usable only if it is complete: the header says how many
%   eigenvalues the run has and how many are written, the two must agree, and the
%   rows must be there. Anything else -- a run interrupted halfway through
%   writing, or a file left over from when only the head was kept -- is a miss
%   rather than a short spectrum, so that a truncated cache can never be read as
%   a run that has no such eigenvalue.
    lambda = NaN; dofs = NaN; tsec = NaN;
    if ~exist(path, 'file')
        status = 'miss';
        return;
    end
    fid = fopen(path, 'r');
    if fid == -1
        status = 'miss';
        return;
    end
    closer = onCleanup(@() fclose(fid));
    stored = NaN; total = NaN;
    while true                      % the header, which is the "#" lines
        line = fgetl(fid);
        if ~ischar(line) || isempty(line) || line(1) ~= '#'; break; end
        tok = regexp(line, 'dofs\s*=\s*(\d+)', 'tokens', 'once');
        if ~isempty(tok); dofs = str2double(tok{1}); end
        tok = regexp(line, 'time:\s*([0-9.eE+-]+)', 'tokens', 'once');
        if ~isempty(tok); tsec = str2double(tok{1}); end
        tok = regexp(line, 'stored:\s*(\d+)\s+of\s+(\d+)', 'tokens', 'once');
        if ~isempty(tok)
            stored = str2double(tok{1});
            total  = str2double(tok{2});
        end
    end
    clear closer;                   % READTABLE opens the file again itself
    % The whole spectrum of a fine run is ten thousand rows, so the table is read
    % in one call rather than line by line, as READ_EIGS_CSV does it.
    try
        evals = readtable(path, 'CommentStyle', '#').lambda_n;
    catch
        status = 'miss';
        return;
    end
    if isnan(stored) || stored ~= total || numel(evals) ~= total
        status = 'miss';       % incomplete, truncated or hand-edited
        return;
    end
    if index <= total
        lambda = evals(index);
        status = 'hit';
    else
        lambda = total;        % the whole spectrum is here and it is too short
        status = 'short';
    end
end


function write_spectrum(path, cfg, p, dofs, tsec, evals)
%WRITE_SPECTRUM One run's whole spectrum, with what it took to get it.
%
%   Every eigenvalue the run resolved, not a head of the spectrum: what the dense
%   eig produced is what is kept, and no index can then be asked for that the
%   cache has to go back to the solver for. It costs a few megabytes over the
%   three domains, against a quarter of an hour of eig per domain.
%
%   Seventeen significant digits, and the header repeats the count the rows carry
%   so that READ_SPECTRUM can tell a complete file from an interrupted one.
    keep = numel(evals);
    fid = fopen(path, 'w');
    if fid == -1
        error('eigenvalues_convergence:spectrum', ...
            'Could not open %s for writing.', path);
    end
    closer = onCleanup(@() fclose(fid));
    fprintf(fid, '# Domain: %s (%s, convergence sweep, dense eig full spectrum)\n', ...
            cfg.name, p.method);
    fprintf(fid, '# Resolution %s, dofs = %d\n', p.resolution, dofs);
    fprintf(fid, '# Computation time: %.4f s\n', tsec);
    fprintf(fid, '# Eigenvalues stored: %d of %d\n', keep, numel(evals));
    fprintf(fid, 'n,lambda_n\n');
    % One call rather than a loop over the rows: the spectra of the finer runs
    % are ten thousand lines long, and there are dozens of them per domain.
    fprintf(fid, '%d,%.17g\n', [1:keep; reshape(evals(1:keep), 1, [])]);
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
    % evals = sort(-real(eig(L)), 'ascend');
    evals = sort(-eig(dst_laplace_symmetrise(L)), 'ascend');
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
%FIND_RUN Whether a run table holds this run.
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
        % eig(L);
        eig(dst_laplace_symmetrise(L));
        fd_laplace_spectrum(struct('box', cfg.box, 'phi', cfg.phi, 'M', M));
        fem_laplace_spectrum(entry, "eig");
    end
end


function s = mk(method, resolution, dofs, lambda, tsec)
    s = struct('method', method, 'resolution', resolution, 'dofs', dofs, ...
               'lambda', lambda, 'time', tsec);
end


function report(r, index, whence)
%REPORT One line per run, saying where it came from when it was not computed here.
    fprintf('  %-4s %-12s dofs = %-5d  lambda_%d = %.8f  time = %6.2f s %s\n', ...
            r.method, r.resolution, r.dofs, index, r.lambda, r.time, whence);
end


function write_runs_csv(csv, cfg, runs, index, ref)
%WRITE_RUNS_CSV One row per run: method, resolution, dofs, eigenvalue, time.
%
%   The reference is named in the header but kept out of the table: it is not a
%   run, and the figures recompute it from the config anyway.
%
%   Seventeen significant digits: the CSV is the cache the figures are redrawn
%   from, so it has to carry a double back exactly. Twelve, which is plenty to
%   read, would round an error that sits near machine precision -- the DST error
%   on the rectangle -- to zero or to a spurious level, and the redrawn figure
%   would not be the one the run produced.
    fid = fopen(csv, 'w');
    if fid == -1
        error('eigenvalues_convergence:csv', 'Could not open %s for writing.', csv);
    end
    closer = onCleanup(@() fclose(fid));
    fprintf(fid, ['# Domain: %s (convergence of eigenvalue %d against DOF, ', ...
                  'dense eig)\n'], cfg.name, index);
    fprintf(fid, '# Reference (%s): lambda_%d = %.13f\n', ref.label, index, ref.value);
    fprintf(fid, 'method,resolution,dofs,lambda_%d,time_s\n', index);
    for i = 1:numel(runs)
        fprintf(fid, '%s,%s,%d,%.17g,%.4f\n', runs(i).method, runs(i).resolution, ...
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


function ref = reference_for(cfg, index, dst_ref)
%REFERENCE_FOR The value the figures measure against, and what it is worth.
%
%   EXACT, where the domain has a closed-form spectrum: the rectangle and the
%   right isosceles triangle, at every index. It is external to the sweep, so
%   every run can be measured against it and the error is a true error.
%
%   PUBLISHED, when the index asked for is the one the config carries a value
%   for: the ground state of the L-shape and of the four catalogue domains,
%   against the MPS values of Betcke & Trefethen. Again external to the sweep,
%   and again a true error.
%
%   Otherwise DST_REF, a DST run of its own at the finest resolution a dense eig
%   fits on this machine, computed by DST_REFERENCE. It is not a true reference
%   -- what it gives is a difference between two computed numbers rather than a
%   distance from the exact value, and being DST itself it shares whatever error
%   the DST runs of the sweep have -- but it lies well beyond the sweep, so the
%   whole of every curve can be drawn against it. That is what it buys over the
%   finest run of the sweep, which the figures used to fall back on and which
%   left the three finest DST points measuring themselves.
%
%   What is left is a floor rather than an error: the finest DST runs approach
%   the reference's own resolution and cannot be told from it. The figures do not
%   mark that floor; the run table names the reference and its dof count.
%
%   Either way the error is reported relative to the reference, so that the
%   figures of different indices can be read against each other.
    if ~isempty(cfg.exact_fun)
        % "exact" reads as a word, not as a method, so it is set upright in the
        % legend where MPS is set in typewriter.
        ref = struct('value', cfg.exact_fun(index), 'label', 'exact', ...
                     'legend', 'exact', ...
                     'published', true, 'relative', true, 'dst_trim', 0);
        fprintf('  reference: exact, lambda_%d = %.13f\n', index, ref.value);
        return;
    end
    if ~isempty(cfg.ref_index) && index == cfg.ref_index
        ref = struct('value', cfg.ref_value, 'label', cfg.ref_label, ...
                     'legend', sprintf('\\texttt{%s}', cfg.ref_label), ...
                     'published', true, 'relative', true, 'dst_trim', 0);
        fprintf('  reference: %s, lambda_%d = %.13f\n', ref.label, index, ref.value);
        return;
    end
    if isnan(dst_ref.value)
        error('eigenvalues_convergence:noReference', ...
            ['No published value for lambda_%d on %s and no reference run: ', ...
             'give the domain an M_ref.'], index, cfg.name);
    end
    % No trimming: the reference is not a member of the sweep.
    ref = struct('value', dst_ref.value, 'label', 'DST', ...
                 'legend', '\texttt{DST}', ...
                 'published', false, 'relative', true, 'dst_trim', 0);
    fprintf('  reference: DST at %d dofs, lambda_%d = %.9f\n', ...
            dst_ref.dofs, index, ref.value);
end


function dst_ref = dst_reference(cfg, index, spectra_dir)
%DST_REFERENCE The DST run the figures measure against where nothing is published.
%
%   One run per domain, at cfg.M_ref, which is the finest grid whose dense eig
%   fits in memory here: about 20000 degrees of freedom, a matrix of 2.9 GB that
%   eig doubles while it works. Past that it swaps, and the first thousand
%   eigenvalues would have to come from eigs on the sparse operator instead --
%   which agrees with the dense answer to 4e-12 and has no memory ceiling, but
%   costs ten minutes at 30000 degrees of freedom against two at 20000.
%
%   It goes through the same spectra cache as the runs of the sweep, so it is
%   computed once per domain however many indices are asked for, and a domain
%   that needs no reference -- one with a closed form, or the ground state
%   against its published value -- never computes it at all.
    dst_ref = struct('value', NaN, 'dofs', NaN);
    if ~isempty(cfg.exact_fun) || (~isempty(cfg.ref_index) && index == cfg.ref_index)
        return;
    end
    if isempty(cfg.M_ref)
        return;
    end

    resolution = sprintf('M = %d', cfg.M_ref);
    path = spectrum_path(spectra_dir, cfg, 'DST', resolution);
    [lambda, dofs, ~, status] = read_spectrum(path, index);
    if strcmp(status, 'hit')
        dst_ref = struct('value', lambda, 'dofs', dofs);
        return;
    end
    if strcmp(status, 'short')
        error('eigenvalues_convergence:referenceShort', ...
            'The reference run of %s has only %d eigenvalues, lambda_%d was asked for.', ...
            cfg.name, lambda, index);
    end

    [x_range, y_range] = bounding_box(cfg.box(1), cfg.box(2), cfg.box(3), cfg.box(4));
    fprintf('  reference run: DST at M = %d, computing ...\n', cfg.M_ref);
    t = tic;
    [evals, dofs] = dst_spectrum(x_range, y_range, cfg.M_ref, cfg.phi);
    tsec = toc(t);
    write_spectrum(path, cfg, struct('method', 'DST', 'resolution', resolution), ...
                   dofs, tsec, evals);
    fprintf('  reference run: %d dofs, %.1f s\n', dofs, tsec);
    if numel(evals) < index
        error('eigenvalues_convergence:referenceShort', ...
            'The reference run of %s has only %d eigenvalues, lambda_%d was asked for.', ...
            cfg.name, numel(evals), index);
    end
    dst_ref = struct('value', evals(index), 'dofs', dofs);
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
            'DisplayName', ref.legend);
    end

    hold(ax, 'off');
    xlabel(ax, 'degrees of freedom');
    ylabel(ax, sprintf('$\\lambda_{%d}$', index));
    % No title: the domain is identified by the output file name.
    % The corner the curves leave free is not the same one from domain to domain
    % -- on the L-shape they come down to the level from above and clear the top
    % right, on the rectangle they come up to it and clear the bottom right -- so
    % the placement is left to the axes rather than fixed here.
    legend(ax, 'Location', 'best', 'FontSize', 10, 'Box', 'off');
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
%   Only the part of a curve above ROUNDOFF_FLOOR is fitted, and a curve with
%   fewer than three points there, or one that does not fall over them (see
%   DECAYS), gets no triangle whatever the placement asks for. Both cases are the
%   DST curve on the rectangle, where the sine basis is exact: at the ground
%   state the whole curve is the eigensolver's roundoff, which grows with the
%   norm of the operator rather than falling with the grid, and deeper in the
%   spectrum it is one coarse grid short of the modes it needs and then drops to
%   that same floor in a single step.
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
    % Where the arithmetic takes over from the method. A dense eig returns the
    % eigenvalues of these operators to about a part in 1e11 -- the norm grows as
    % h^-2, and that is what the DST error on the rectangle is made of -- so a
    % relative error a decade below that is not a discretisation error and no
    % rate is fitted through it. Every genuine curve in these figures is orders
    % of magnitude above it, the closest being FEM on the rectangle ground state
    % at 1e-7.
    ROUNDOFF_FLOOR = 1e-10;
    % How many of the finest runs the rate is fitted to; see the fit below.
    N_FIT = 4;

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
        % The spans are fractions of the drawing window, the finer half of the
        % curve, so all three triangles live at the end the rates belong to;
        % they are staggered across it so that no two stack up. FEM keeps the shortest of them: it
        % is the lowest curve and the steepest, and a taller triangle would grow
        % past the gap to the curve above it.
        spans    = {[0.38 0.66], [0.05 0.33], [0.70 0.98]};
        offsets  = [1.13, 1.25, 1.10];
    else
        triangle = [false, true, true];
        spans    = {[], [0.05 0.45], [0.55 0.95]};
        offsets  = [NaN, 2.2, 1.8];
    end

    lo = inf; hi = 0;
    for mi = 1:numel(methods)
        [d, v] = curve(runs, methods{mi}, ref.dst_trim);
        err = abs(v - ref.value);
        if ref.relative
            err = err / abs(ref.value);
        end
        % A run that lands on the reference exactly -- it happens on the
        % rectangle, where DST is exact -- has no place on a log axis, and
        % dropping it would leave a hole in the middle of a curve. It is drawn on
        % the floor instead, at the smallest error double precision can tell from
        % zero.
        err = max(err, eps);
        loglog(ax, d, err, ['-' markers{mi}], 'Color', colors{mi}, 'LineWidth', 2.0, ...
            'MarkerSize', 7, 'MarkerFaceColor', 'w', ...
            'DisplayName', sprintf('\\texttt{%s}', methods{mi}));
        yhi = 0;
        % The rate is fitted to the finest N_FIT runs of the part of the curve
        % that is a discretisation error, which is the part above the floor.
        %
        % Above the floor, because on the rectangle DST at an index deep enough
        % to feel the edge of the coarsest grid drops from a real error to
        % roundoff in one step, and a line through the drop describes neither end
        % of it.
        %
        % The finest runs, because a rate is an asymptotic statement and the
        % coarse end of a curve is not asymptotic. FD on the GWW drums is the
        % case that forces it: its coarsest run happens to land within 4e-4 of
        % the ground state, the curve then rises before it falls, and a line
        % through the whole of it reported -0.09 -- a description of the hump and
        % not of the method. The count follows the DOF-requirement family, which
        % fits its extrapolations the same way.
        % Two questions, two windows. Whether the curve converges at all is asked
        % of the whole of it, since the finest few runs of a slow method fall by
        % a fifth and would fail a test meant to catch curves that do not fall.
        % At what rate is asked of the finest runs alone.
        above = find(err > ROUNDOFF_FLOOR);
        fit   = above(max(1, numel(above) - N_FIT + 1) : end);
        if triangle(mi) && numel(above) >= 3 && numel(fit) >= 3 && decays(err(above))
            % Least-squares algebraic rate: err ~ dofs^p, annotated by the
            % triangle, which is drawn over the runs the rate was fitted to.
            p = polyfit(log(d(fit)), log(err(fit)), 1);
            % Drawn over the finer half of the curve rather than over the four
            % runs the rate was fitted to: four runs of a sweep span a fifth of a
            % decade, and three triangles crammed into that are unreadable. The
            % half is still the end the rate belongs to.
            draw = above(ceil(numel(above) / 2) : end);
            yhi = slope_triangle(ax, d(draw), err(draw), p, colors{mi}, ...
                                 spans{mi}, offsets(mi));
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
    % No title: the domain is identified by the output file name. The free corner
    % moves with the domain here too, so the placement is left to the axes; see
    % PLOT_RUNS.
    legend(ax, 'Location', 'best', 'FontSize', 10, 'Box', 'off');
    grid(ax, 'on'); box(ax, 'on');
    set(ax, 'FontSize', 12);

    save_figure(fig, stem);
end


function tf = decays(err)
%DECAYS Whether an error curve falls fast enough for a fitted rate to mean anything.
%
%   The finest run at least twice as accurate as the worst, over at least three
%   runs. Any method converging at all clears that -- the shallowest curve in
%   these figures, FD on the L-shaped ground state, falls by a factor of eight --
%   while a curve that wanders about one level or climbs does not: the DST error
%   on the rectangle is the eigensolver's roundoff, which grows with the norm of
%   the operator, and a line fitted to it would report a rate for arithmetic
%   noise. Written as a ratio rather than as a threshold on the error itself, so
%   that it carries no scale of its own.
    tf = numel(err) >= 3 && err(end) < 0.5 * max(err);
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


function y_hi = slope_triangle(ax, d, err, p, color, span, offset)
%SLOPE_TRIANGLE Reference-slope triangle above one curve, labelled with the rate.
%
%   Right triangle whose hypotenuse has the fitted slope p(1), drawn just above
%   the curve it belongs to: horizontal leg on top, vertical leg on the right.
%   SPAN gives the fraction of the log-DOF range it covers and OFFSET how far
%   above the curve it sits -- the caller staggers the three triangles so that
%   each keeps to the gap above its own curve without running into the
%   neighbouring one. Returns the highest ordinate drawn so the caller can
%   leave room. Kept out of the legend.
%
%   Anchored to the curve at the left end of the span, not to the fitted line.
%   The two agree wherever the curve is a power law, which is most of these
%   figures; where it is not -- DST on the triangle, steep while the coarse grids
%   are still resolving the mode and shallower afterwards -- the fitted line runs
%   away from the curve in the middle and a triangle hung off it floats in empty
%   space. The slope drawn is the fitted one either way.

    lg = log(d([1 end]));
    x1 = exp(lg(1) + span(1) * diff(lg));
    x2 = exp(lg(1) + span(2) * diff(lg));
    y1 = exp(interp1(log(d), log(err), log(x1))) * offset;
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
