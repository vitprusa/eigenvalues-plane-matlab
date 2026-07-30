function make_eigenvalues_head_paper_tables(name)
%MAKE_EIGENVALUES_HEAD_PAPER_TABLES Paper eigenvalue-comparison tables (self-computed).
%
%   Paper (article) counterpart of
%   experiments/eigenvalues_head/make_eigenvalues_head_comparison_tables.py.
%   Unlike that script -- which reads the committed per-method CSVs -- this one
%   GENERATES ITS OWN DATA: for each domain it recomputes the Dirichlet-Laplacian
%   spectrum with DST, finite differences, finite elements (all dense eig) and,
%   for the rectangular domains, Chebyshev spectral collocation, at the first
%   few resolutions of each method, timing every run with tic/toc.
%
%   For each domain it writes
%     - one CSV per (method, resolution) into results_paper/eigenvalues_head/
%       (a NEW location; the existing results/eigenvalues[_dof_sweep]/ CSVs are
%       never touched), each carrying the dof count and the measured time in its
%       metadata header, and
%     - one LaTeX snippet <domain>_eigenvalues_head.tex holding a booktabs
%       table: rows are the ground-truth reference plus each method's DOF runs,
%       columns are DOF, Time (s) and the first eight eigenvalues, and
%     - <domain>_eigenvalues_head_transposed.tex, the same data with the axes
%       swapped -- one column per run, one row per eigenvalue index, so that a
%       single eigenvalue reads across all discretisations, and
%     - <domain>_eigenvalues_head_transposed_extended.tex for the domains whose
%       config carries EXTRA_INDICES: the transposed table with a few deeper
%       eigenvalues appended, each on an empty row of its own (see
%       WRITE_LATEX_TRANSPOSED).
%   All are \scriptsize and print four decimals, so that the tables of a domain
%   match; the extended one is too wide even so, and is scaled to the text width
%   with \resizebox.
%
%   The ground truth is the analytic spectrum where a closed form is known
%   (rectangle, isosceles triangle) or the MPS reference (all other domains);
%   it is read as a reference, not recomputed. The MPS references of the
%   non-analytic domains come from data/, which holds only the leading few
%   eigenvalues -- the remaining reference cells are left blank.
%
%   Scope is seven domains, in two groups. The three DOF-sweep domains
%   (rectangle, isosceles_triangle, L_shaped) -- those for which a four-level,
%   cross-method DOF schedule is defined (see
%   experiments/eigenvalues_dof_sweep/compute_*_dof_sweep.m) -- run at four DOF
%   levels. The four non-analytic domains (ellipse_minus_quadrant, H_shaped,
%   gww1, gww2) have no DOF sweep of their own; their levels are set locally in
%   DOMAIN_CONFIGS -- four grid resolutions for DST/FD, three meshes for FEM.
%
%   Called with no argument (or an empty one) it does every domain; pass a
%   domain name to regenerate just that one.
%
%   Requires the PDE Toolbox (FEM) and Chebfun (Cheb).

    NEIG = 8;    % eigenvalues shown per table

    here         = fileparts(mfilename('fullpath'));
    project_root = fileparts(fileparts(here));
    run(fullfile(project_root, 'startup.m'));
    addpath(fullfile(project_root, 'experiments'));   % domain_catalog_dst / _fem

    out_dir = fullfile(project_root, 'results_paper', 'eigenvalues_head');
    if ~exist(out_dir, 'dir')
        mkdir(out_dir);
    end

    cfgs = domain_configs(project_root);
    if nargin >= 1 && ~isempty(name)
        cfgs = cfgs(strcmp({cfgs.name}, name));
        if isempty(cfgs)
            error('head_paper:unknownDomain', ...
                ['Unknown domain "%s" (known: rectangle, isosceles_triangle, ', ...
                 'L_shaped, ellipse_minus_quadrant, H_shaped, gww1, gww2).'], name);
        end
    end
    for i = 1:numel(cfgs)
        cfg = cfgs(i);
        fprintf('=== %s ===\n', cfg.name);
        rows = compute_domain(cfg, out_dir);
        tex  = fullfile(out_dir, sprintf('%s_eigenvalues_head.tex', cfg.name));
        write_latex(tex, cfg, rows, NEIG);
        fprintf('Wrote %s\n', tex);
        tex = fullfile(out_dir, sprintf('%s_eigenvalues_head_transposed.tex', cfg.name));
        write_latex_transposed(tex, cfg, rows, NEIG, []);
        fprintf('Wrote %s\n', tex);
        if ~isempty(cfg.extra_indices)
            tex = fullfile(out_dir, ...
                sprintf('%s_eigenvalues_head_transposed_extended.tex', cfg.name));
            write_latex_transposed(tex, cfg, rows, NEIG, cfg.extra_indices);
            fprintf('Wrote %s\n', tex);
        end
    end
end


function cfgs = domain_configs(project_root)
%DOMAIN_CONFIGS Per-domain geometry, DOF schedules and ground-truth source.
%
%   EXTRA_INDICES holds the deeper eigenvalue indices of the extended transposed
%   table, empty for the domains that do not get one.
    mps = fullfile(project_root, 'results', 'eigenvalues', 'mps', ...
                   'L_shaped_eigenvalues_MPS.csv');

    cfgs = struct('name', {}, 'pretty', {}, 'box', {}, 'phi', {}, ...
                  'gd', {}, 'ns', {}, 'sf', {}, 'M_grid', {}, 'Hmax_fem', {}, ...
                  'has_cheb', {}, 'N_cheb', {}, 'truth_kind', {}, ...
                  'analytic_fun', {}, 'mps_path', {}, 'extra_indices', {});

    % --- rectangle [0, 2*pi] x [0, pi] --------------------------------------
    % The only domain with an extended table so far: its ground truth is a
    % closed form, so the deeper eigenvalues have an exact value to be read
    % against, which the short MPS references of the other domains cannot give.
    cfgs(end+1) = struct( ...
        'name', 'rectangle', 'pretty', 'rectangular', ...
        'box', [0 2*pi 0 pi], ...
        'phi', @(x, y) indicator_rectangle(x, y, 0, 2*pi, 0, pi), ...
        'gd', [3; 4; 0; 2*pi; 2*pi; 0; 0; 0; pi; pi], 'ns', char('R1')', 'sf', 'R1', ...
        'M_grid', [15 25 35 49], 'Hmax_fem', [0.35 0.25 0.18 0.13], ...
        'has_cheb', true, 'N_cheb', [10 20 30 40], ...
        'truth_kind', 'analytic', ...
        'analytic_fun', @() analytic_rectangle(), 'mps_path', '', ...
        'extra_indices', [60 80 100 200]);

    % --- right isosceles triangle, legs pi ----------------------------------
    cfgs(end+1) = struct( ...
        'name', 'isosceles_triangle', 'pretty', 'right isosceles triangle', ...
        'box', [0 pi 0 pi], ...
        'phi', @(x, y) indicator_isosceles_triangle(x, y, 0, pi, 0), ...
        'gd', [2; 3; 0; pi; pi; 0; 0; pi], 'ns', char('T1')', 'sf', 'T1', ...
        'M_grid', [15 25 35 50], 'Hmax_fem', [0.20 0.13 0.09 0.06], ...
        'has_cheb', false, 'N_cheb', [], ...
        'truth_kind', 'analytic', ...
        'analytic_fun', @() analytic_isosceles(), 'mps_path', '', ...
        'extra_indices', []);

    % --- L-shaped domain ----------------------------------------------------
    cfgs(end+1) = struct( ...
        'name', 'L_shaped', 'pretty', 'L-shaped', ...
        'box', [-1 1 -1 1], ...
        'phi', @(x, y) indicator_L_shaped(x, y, -1, 0, 1, -1, 0, 1), ...
        'gd', [[3; 4; -1; 1; 1; -1; -1; -1; 1; 1], [3; 4; 0; 1; 1; 0; 0; 0; 1; 1]], ...
        'ns', char('R1', 'R2')', 'sf', 'R1-R2', ...
        'M_grid', [15 25 35 49], 'Hmax_fem', [0.20 0.13 0.09 0.06], ...
        'has_cheb', false, 'N_cheb', [], ...
        'truth_kind', 'mps', 'analytic_fun', [], 'mps_path', mps, ...
        'extra_indices', []);

    % --- Non-analytic domains: ellipse-minus-quadrant, H, GWW1, GWW2 ---------
    % No closed-form spectrum and no Cheb (non-rectangular), so the table
    % compares DST/FD at four grid levels and FEM at three mesh levels against
    % the MPS reference eigenvalues of Betcke & Trefethen kept in data/ -- only
    % the leading few eigenvalues are published there, so the reference row is
    % short. Box and indicator come from DOMAIN_CATALOG_DST, the FEM decsg
    % geometry from DOMAIN_CATALOG_FEM; the DST/FD grid resolutions M keep the
    % domain edges on grid lines (M+1 divisible by 4 for the ellipse, by 3 for
    % H, by 6 for the 6-wide GWW box). The fourth level is sized to land near
    % 3000 dofs, comparable to the finest FEM mesh.
    % Columns: catalog name (for the DST/FEM lookups), output name (file names,
    % \texttt label and \label -- H uses H_shaped), pretty caption name,
    % DST/FD grid resolutions, FEM mesh sizes, MPS reference CSV in data/.
    extra = { ...
        'ellipse_minus_quadrant', 'ellipse_minus_quadrant', 'ellipse-minus-quadrant', [23 35 47 103], [0.20 0.13 0.09], 'ellipse_minus_quadrant.csv'; ...
        'H',                      'H_shaped',               'H-shaped',              [20 35 50 62],  [0.18 0.12 0.08], 'H_shaped.csv'; ...
        'gww1',                   'gww1',                   'GWW1 isospectral drum', [23 35 47 89],  [0.30 0.20 0.13], 'gww1.csv'; ...
        'gww2',                   'gww2',                   'GWW2 isospectral drum', [23 35 47 89],  [0.30 0.20 0.13], 'gww2.csv'  ...
    };
    for i = 1:size(extra, 1)
        cat_nm = extra{i, 1};
        out_nm = extra{i, 2};
        dstc = domain_catalog_dst(cat_nm);
        femc = domain_catalog_fem(cat_nm);
        cfgs(end+1) = struct( ...
            'name', out_nm, 'pretty', extra{i, 3}, ...
            'box', dstc.box, 'phi', dstc.phi, ...
            'gd', femc.gd, 'ns', femc.ns, 'sf', femc.sf, ...
            'M_grid', extra{i, 4}, 'Hmax_fem', extra{i, 5}, ...
            'has_cheb', false, 'N_cheb', [], ...
            'truth_kind', 'mps', 'analytic_fun', [], ...
            'mps_path', fullfile(project_root, 'data', extra{i, 6}), ...
            'extra_indices', []); %#ok<AGROW>
    end
end


function rows = compute_domain(cfg, out_dir)
%COMPUTE_DOMAIN Compute every method/run for one domain; write CSVs; collect rows.
    [x_range, y_range] = bounding_box(cfg.box(1), cfg.box(2), cfg.box(3), cfg.box(4));

    rows = struct('method_label', {}, 'group', {}, 'dof', {}, 'time', {}, ...
                  'eigs', {}, 'bold', {});

    % Ground-truth reference row (analytic closed form or MPS reference).
    %
    % Each row keeps the whole spectrum, not just the NEIG leading values: the
    % extended table reads deeper indices out of the same rows. The writers take
    % what they need and mark anything past the end of a row with "---", which is
    % how a run coarser than the index asked for, and a published MPS reference
    % shorter than NEIG, are both handled.
    switch cfg.truth_kind
        case 'analytic'
            rows(end+1) = mk('Analytic', 'truth', '---', '---', ...
                             cfg.analytic_fun(), true);
        case 'mps'
            rows(end+1) = mk('\texttt{MPS}', 'truth', '---', '---', ...
                             read_two_col(cfg.mps_path), true);
    end

    methods = {'DST', 'dst'; 'FD', 'fd'; 'FEM', 'fem'};
    if cfg.has_cheb
        methods(end+1, :) = {'Cheb', 'cheb'};
    end

    for mi = 1:size(methods, 1)
        disp_m = methods{mi, 1};
        low    = methods{mi, 2};
        % DST and FD share the grid schedule M_grid; FEM and Cheb walk their
        % own, which need not have the same number of levels.
        switch low
            case 'fem';  nlev = numel(cfg.Hmax_fem);
            case 'cheb'; nlev = numel(cfg.N_cheb);
            otherwise;   nlev = numel(cfg.M_grid);
        end
        for k = 1:nlev
            csv = fullfile(out_dir, sprintf('%s_%s_%d-eigenvalues.csv', ...
                           cfg.name, low, k));
            if exist(csv, 'file')
                % Reuse existing data (e.g. layout-only regeneration): read the
                % eigenvalues, dof count and measured time back, no recompute.
                [evals, dofs, tsec] = read_head_csv(csv);
                fprintf('  %-5s DOF#%d  (reusing %s)\n', disp_m, k, csv);
            else
                try
                    [evals, dofs, resstr, tsec] = run_one(low, k, cfg, x_range, y_range);
                catch ME
                    warning('head_paper:run', '%s %s DOF#%d failed: %s', ...
                            cfg.name, disp_m, k, ME.message);
                    continue;
                end
                write_head_csv(csv, cfg.name, disp_m, resstr, dofs, tsec, evals);
                fprintf('  %-5s DOF#%d  dofs = %-5d  time = %6.2f s  ->  %s\n', ...
                        disp_m, k, dofs, tsec, csv);
            end
            rows(end+1) = mk(disp_m, disp_m, num2str(dofs), sprintf('%.2f', tsec), ...
                             evals, false); %#ok<AGROW>
        end
    end
end


function [evals, dofs, resstr, tsec] = run_one(low, k, cfg, x_range, y_range)
%RUN_ONE One timed spectrum computation for method `low` at DOF level k.
%
%   The timing covers COMPUTE_ONE only, and never a cold call: warm-up is spent
%   on discarded runs by ENSURE_WARM instead. Left in, it lands on whichever run
%   happens to go first and makes a coarse grid look slower than a finer one.
    ensure_warm(low, cfg, x_range, y_range);
    t = tic;
    [evals, dofs, resstr] = compute_one(low, k, cfg, x_range, y_range);
    tsec = toc(t);
end


function ensure_warm(low, cfg, x_range, y_range)
%ENSURE_WARM Discarded runs of method `low` on domain `cfg`, before its first timing.
%
%   Warm-up has two components, both measured on this code: MATLAB's toolbox
%   load and JIT compilation, paid once per session (about 1.7 s for FEM), and a
%   smaller per-geometry cost, paid again by the first call on each new domain
%   (about 0.14 s against 0.03 s warm, for the L-shape). Hence the key is method
%   AND domain, and two calls rather than one -- the first geometry of a session
%   is still decaying on its second call.
%
%   Lazy on purpose: a run whose CSV is cached never reaches RUN_ONE, so a
%   layout-only regeneration still costs no computation at all. The discarded
%   runs use level 1, the cheapest of each sweep.
    persistent warmed
    if isempty(warmed)
        warmed = struct();
    end
    key = sprintf('%s_%s', low, cfg.name);
    if isfield(warmed, key)
        return;
    end
    warmed.(key) = true;   % set first, so that a failed warm-up is not retried
    fprintf('  warm-up: %s on %s (discarded)\n', low, cfg.name);
    try
        for i = 1:2
            compute_one(low, 1, cfg, x_range, y_range);
        end
    catch
        % Ignored: the real run reports its own failure through the caller.
    end
end


function [evals, dofs, resstr] = compute_one(low, k, cfg, x_range, y_range)
%COMPUTE_ONE One untimed spectrum computation for method `low` at DOF level k.
    switch low
        case 'dst'
            M = cfg.M_grid(k);
            resstr = sprintf('M = %d', M);
            [L, info] = make_dst_laplace_mat_batched(x_range, y_range, M, cfg.phi);
            evals = sort(-real(eig(L)), 'ascend');
            dofs = info.dofs;
        case 'fd'
            M = cfg.M_grid(k);
            resstr = sprintf('M = %d', M);
            [evals, info] = fd_laplace_spectrum(struct('box', cfg.box, 'phi', cfg.phi, 'M', M));
            dofs = info.dofs;
        case 'fem'
            h = cfg.Hmax_fem(k);
            resstr = sprintf('Hmax = %g', h);
            entry = struct('gd', cfg.gd, 'ns', cfg.ns, 'sf', cfg.sf, ...
                           'Hmax_eig', h, 'Hmax_solvepdeeig', h);
            [evals, info] = fem_laplace_spectrum(entry, "eig");
            dofs = info.dofs;
        case 'cheb'
            N = cfg.N_cheb(k);
            resstr = sprintf('N = %d', N);
            [evals, info] = chebfun_laplace_spectrum(cfg.box, N);
            dofs = info.dofs;
        otherwise
            error('head_paper:method', 'unknown method %s', low);
    end
    evals = sort(real(evals(:)), 'ascend');
end


function s = mk(method_label, group, dof, time, eigs, bold)
    s = struct('method_label', method_label, 'group', group, 'dof', dof, ...
               'time', time, 'eigs', eigs(:)', 'bold', bold);
end


function [evals, dofs, tsec] = read_head_csv(path)
%READ_HEAD_CSV Read a paper-head CSV back: eigenvalues, dof count, measured time.
    fid = fopen(path, 'r');
    if fid == -1
        error('head_paper:csvread', 'Could not open %s.', path);
    end
    closer = onCleanup(@() fclose(fid));
    dofs = NaN; tsec = NaN; evals = [];
    while true
        line = fgetl(fid);
        if ~ischar(line); break; end
        s = strtrim(line);
        if isempty(s); continue; end
        if s(1) == '#'
            m = regexp(s, 'dofs\s*=\s*(\d+)', 'tokens', 'once');
            if ~isempty(m); dofs = str2double(m{1}); end
            m = regexp(s, 'Computation time:\s*([0-9.]+)', 'tokens', 'once');
            if ~isempty(m); tsec = str2double(m{1}); end
            continue;
        end
        if strncmpi(s, 'n,', 2); continue; end
        parts = strsplit(s, ',');
        if numel(parts) >= 2
            v = str2double(parts{2});
            if ~isnan(v); evals(end+1, 1) = v; end %#ok<AGROW>
        end
    end
end


function v = analytic_rectangle()
%ANALYTIC_RECTANGLE lambda = m^2/4 + n^2 on [0,2pi]x[0,pi], ascending.
    [Mm, Nn] = meshgrid(1:60, 1:40);
    v = sort(reshape(Mm.^2 / 4 + Nn.^2, [], 1), 'ascend');
end


function v = analytic_isosceles()
%ANALYTIC_ISOSCELES lambda = m^2 + n^2 with m > n >= 1 (legs pi), ascending.
    [Mm, Nn] = meshgrid(1:60, 1:60);
    keep = Mm > Nn;
    v = sort(Mm(keep).^2 + Nn(keep).^2, 'ascend');
end


function v = read_two_col(path)
%READ_TWO_COL Second column of a "# ... / n,lambda_n / i,value" CSV.
    fid = fopen(path, 'r');
    if fid == -1
        error('head_paper:mps', 'Could not open reference CSV %s.', path);
    end
    closer = onCleanup(@() fclose(fid));
    v = [];
    while true
        line = fgetl(fid);
        if ~ischar(line); break; end
        s = strtrim(line);
        if isempty(s) || s(1) == '#'; continue; end
        parts = strsplit(s, ',');
        if numel(parts) < 2; continue; end
        val = str2double(parts{2});
        if ~isnan(val)
            v(end+1, 1) = val; %#ok<AGROW>
        end
    end
end


function write_head_csv(csv, name, method, resstr, dofs, tsec, evals)
%WRITE_HEAD_CSV One spectrum CSV with dof/time metadata, then n,lambda_n rows.
    fid = fopen(csv, 'w');
    if fid == -1
        error('head_paper:csv', 'Could not open %s for writing.', csv);
    end
    fprintf(fid, '# Domain: %s (%s, paper eigenvalues-head run, dense eig full spectrum)\n', ...
            name, method);
    fprintf(fid, '# Resolution %s, dofs = %d\n', resstr, dofs);
    fprintf(fid, '# Computation time: %.4f s\n', tsec);
    fprintf(fid, 'n,lambda_n\n');
    for i = 1:numel(evals)
        fprintf(fid, '%d,%.12g\n', i, evals(i));
    end
    fclose(fid);
end


function write_latex(tex, cfg, rows, NEIG)
%WRITE_LATEX Transposed booktabs table for one domain.
    fid = fopen(tex, 'w');
    if fid == -1
        error('head_paper:tex', 'Could not open %s for writing.', tex);
    end
    closer = onCleanup(@() fclose(fid));

    t = caption_bits(cfg);
    cheb_note         = t.cheb_note;
    leading           = t.leading;
    comparison_clause = t.comparison_clause;
    lev_note          = t.lev_note;
    method_list       = t.method_list;

    colspec = ['lrr', repmat('r', 1, NEIG)];

    % Header comment.
    fprintf(fid, ['%% Eigenvalue comparison for the %s domain \\texttt{%s}: DST, FD, ', ...
        'FEM%s,\n'], cfg.pretty, latex_name(cfg.name), cheb_note);
    fprintf(fid, '%% %s, first %d eigenvalues, self-computed with timing.\n', lev_note, NEIG);
    fprintf(fid, '%%\n%% Requires in the preamble:\n');
    fprintf(fid, ['%%   \\usepackage{booktabs}\n%%   \\usepackage{multirow}\n', ...
                  '%%   \\usepackage{amsmath}\n%%   \\usepackage{listings}\n\n']);

    fprintf(fid, '\\begin{table}[htbp]\n  \\centering\n');
    % \scriptsize scopes to the tabular (grouped) so the caption stays at normal
    % size; both layouts use it so that the two tables of a domain match.
    fprintf(fid, '  {\\scriptsize\n  \\begin{tabular}{%s}\n', colspec);
    fprintf(fid, '    \\toprule\n');
    fprintf(fid, '    Method & DOF & Time (s)');
    for i = 1:NEIG
        fprintf(fid, ' & $\\lambda_{%d}$', i);
    end
    fprintf(fid, ' \\\\\n    \\midrule\n');

    % Rows are grouped by method (contiguous). The method name is written once
    % per group with \multirow spanning the group's rows; the individual DOF
    % runs are distinguished by the DOF column.
    r = 1;
    while r <= numel(rows)
        grp = rows(r).group;
        j = r;
        while j <= numel(rows) && strcmp(rows(j).group, grp)
            j = j + 1;
        end
        nrows = j - r;
        if r > 1
            fprintf(fid, '    \\midrule\n');
        end
        for idx = r:(j - 1)
            row = rows(idx);
            if idx == r
                if strcmp(grp, 'truth')
                    firstcell = row.method_label;
                else
                    firstcell = sprintf('\\multirow{%d}{*}{\\texttt{%s}}', nrows, row.method_label);
                end
            else
                firstcell = '';
            end
            fprintf(fid, '    %s & %s & %s', firstcell, row.dof, row.time);
            for i = 1:NEIG
                if i <= numel(row.eigs)
                    s = sprintf('%.4f', row.eigs(i));
                    if row.bold
                        s = sprintf('\\textbf{%s}', s);
                    end
                else
                    % Not available (e.g. the published MPS references list
                    % fewer than NEIG eigenvalues); same marker as the DOF and
                    % Time cells of the reference row.
                    s = '---';
                end
                fprintf(fid, ' & %s', s);
            end
            fprintf(fid, ' \\\\\n');
        end
        r = j;
    end

    fprintf(fid, '    \\bottomrule\n  \\end{tabular}}\n');
    fprintf(fid, ['  \\caption{%s eigenvalues of Laplace operator on the %s ', ...
        'domain \\texttt{%s}. Numerical eigenvalues computed by each ', ...
        'discretisation %s with varying degrees of freedom (DOF)%s. Timing ', ...
        'shows the time for the matrix assembly and full spectrum computation ', ...
        'using \\texttt{MATLAB}''s \\lstinline{eig} function with the default ', ...
        'settings.}\n'], ...
        leading, cfg.pretty, latex_name(cfg.name), method_list, comparison_clause);
    fprintf(fid, '  \\label{tab:eigenvalues_head_%s}\n', cfg.name);
    fprintf(fid, '\\end{table}\n');
end


function t = caption_bits(cfg)
%CAPTION_BITS Shared caption/header wording for both table layouts.
    if cfg.has_cheb
        t.cheb_note = ', Cheb';
    else
        t.cheb_note = '';
    end
    switch cfg.truth_kind
        case 'analytic'
            t.leading = 'Analytical/numerical';
            t.comparison_clause = '; compared with the analytic eigenvalues (exact)';
        case 'mps'
            t.leading = 'Analytical/numerical';
            t.comparison_clause = ['; compared with the reference eigenvalues from the ', ...
                'method of particular solutions (\texttt{MPS})'];
        otherwise   % 'none': no closed-form or reference spectrum available
            t.leading = 'Numerical';
            t.comparison_clause = '';
    end
    nlev = numel(cfg.M_grid);
    nfem = numel(cfg.Hmax_fem);
    if nfem == nlev
        t.lev_note = sprintf('%d DOF runs each', nlev);
    else
        t.lev_note = sprintf('%d DOF runs each (FEM: %d)', nlev, nfem);
    end

    % Method list in \texttt, e.g. "\texttt{DST}, \texttt{FD}, \texttt{FEM} and
    % \texttt{Cheb}" (Cheb only for the rectangular domains).
    mnames = {'DST', 'FD', 'FEM'};
    if cfg.has_cheb
        mnames{end+1} = 'Cheb';
    end
    tt = cellfun(@(m) sprintf('\\texttt{%s}', m), mnames, 'UniformOutput', false);
    if numel(tt) == 1
        t.method_list = tt{1};
    else
        t.method_list = sprintf('%s and %s', strjoin(tt(1:end-1), ', '), tt{end});
    end
end


function write_latex_transposed(tex, cfg, rows, NEIG, extras)
%WRITE_LATEX_TRANSPOSED Transposed booktabs table: one row per eigenvalue.
%
%   The same data as WRITE_LATEX with the axes swapped -- one column per run
%   (the reference plus every method/DOF level), one row per eigenvalue index,
%   so that a single eigenvalue can be read across all discretisations. The
%   method name spans its DOF runs as a \multicolumn group head, the DOF and
%   timing become the two leading body rows, and the reference column stays
%   bold. Four decimals and \scriptsize, as in WRITE_LATEX, which is also what
%   keeps the thirteen to seventeen numeric columns inside the text width.
%
%   EXTRAS are eigenvalue indices past the leading NEIG, each written on a row
%   of its own after an empty row, so that the jump in the index is visible. Pass
%   [] for the plain table; a non-empty EXTRAS also wraps the tabular in
%   \resizebox, drops the timing sentence from the caption and marks the \label
%   as extended, every layout being \input into one document.
    fid = fopen(tex, 'w');
    if fid == -1
        error('head_paper:tex', 'Could not open %s for writing.', tex);
    end
    closer = onCleanup(@() fclose(fid));

    t = caption_bits(cfg);
    ncol = numel(rows);
    colspec = ['l', repmat('r', 1, ncol)];

    % Contiguous column groups: the single truth column, then one per method.
    starts = 1;
    for i = 2:ncol
        if ~strcmp(rows(i).group, rows(i - 1).group)
            starts(end + 1) = i; %#ok<AGROW>
        end
    end
    stops = [starts(2:end) - 1, ncol];

    % Header comment.
    fprintf(fid, ['%% Eigenvalue comparison for the %s domain \\texttt{%s}: DST, FD, ', ...
        'FEM%s,\n'], cfg.pretty, latex_name(cfg.name), t.cheb_note);
    fprintf(fid, ['%% %s, first %d eigenvalues, self-computed with timing.\n', ...
        '%% Transposed layout: one column per run, one row per eigenvalue.\n'], ...
        t.lev_note, NEIG);
    if ~isempty(extras)
        fprintf(fid, ['%% Extended with the eigenvalues %s;\n', ...
            '%% too wide for the text width, so the tabular is scaled to it.\n'], ...
            index_list(extras));
    end
    fprintf(fid, '%%\n%% Requires in the preamble:\n');
    fprintf(fid, ['%%   \\usepackage{booktabs}\n', ...
                  '%%   \\usepackage{amsmath}\n%%   \\usepackage{listings}\n']);
    if ~isempty(extras)
        fprintf(fid, '%%   \\usepackage{graphicx}          %% \\resizebox\n');
    end
    fprintf(fid, '\n');

    % \scriptsize throughout, as in WRITE_LATEX. Measured with pdflatex against
    % the article preamble (amsart, a4paper, geometry scale=0.9, so \textwidth =
    % 537.75pt): at that size 3pt column separation is the widest that still
    % holds the seventeen runs of the rectangle -- the only domain with Cheb --
    % inside the text width, so every table takes it and they stay uniform. A
    % narrower page (a plain article \textwidth of 345pt, say) fits none of them.
    %
    % The extended table is the one that overruns, its deeper eigenvalues running
    % into three figures before the decimal point: a column being as wide as its
    % widest cell, those rows alone set every column width. Four decimals
    % throughout is deliberate -- the table is not to print a different number of
    % digits for the deeper eigenvalues -- and at four decimals no setting of
    % size and separation fits: 105pt over at \scriptsize/3pt, 63pt at \tiny/3pt,
    % 33pt at \scriptsize/1pt. It therefore keeps the family's \scriptsize and
    % 3pt, like every other table here, and buys the fit by scaling the whole
    % tabular to \textwidth with \resizebox (about 84 per cent).
    fprintf(fid, '\\begin{table}[htbp]\n  \\centering\n');
    % The size and \tabcolsep changes scope to the tabular (grouped) so that the
    % caption keeps the normal body size.
    fprintf(fid, '  {\\scriptsize\n  \\setlength{\\tabcolsep}{3pt}\n');
    if ~isempty(extras)
        % Trailing %% so that the line break adds no space before the tabular.
        fprintf(fid, '  \\resizebox{\\textwidth}{!}{%%\n');
    end
    fprintf(fid, '  \\begin{tabular}{%s}\n    \\toprule\n', colspec);

    % Group head: method name spanning its DOF runs, underlined by \cmidrule.
    fprintf(fid, '   ');
    for g = 1:numel(starts)
        lab = rows(starts(g)).method_label;
        if ~strcmp(rows(starts(g)).group, 'truth')
            lab = sprintf('\\texttt{%s}', lab);
        end
        fprintf(fid, ' & \\multicolumn{%d}{c}{%s}', stops(g) - starts(g) + 1, lab);
    end
    fprintf(fid, ' \\\\\n   ');
    for g = 1:numel(starts)
        fprintf(fid, ' \\cmidrule(lr){%d-%d}', starts(g) + 1, stops(g) + 1);
    end
    fprintf(fid, '\n');

    % Run metadata, then the eigenvalues.
    fprintf(fid, '    DOF');
    fprintf(fid, ' & %s', rows.dof);
    fprintf(fid, ' \\\\\n    Time (s)');
    fprintf(fid, ' & %s', rows.time);
    fprintf(fid, ' \\\\\n    \\midrule\n');

    for i = 1:NEIG
        write_eig_row(fid, rows, i, 4);
    end
    % The deeper eigenvalues, each after an empty row: the index jumps from one
    % to the next, and a rule would read as a new block of the same sequence.
    for i = extras(:)'
        fprintf(fid, '   %s \\\\\n', repmat(' &', 1, ncol));
        write_eig_row(fid, rows, i, 4);
    end

    if isempty(extras)
        fprintf(fid, '    \\bottomrule\n  \\end{tabular}}\n');
    else
        fprintf(fid, '    \\bottomrule\n  \\end{tabular}}}\n');   % tabular, \resizebox, size group
    end
    % The extended table stops after the comparison clause: the timing sentence
    % is carried by the two tables it shares its numbers with, and the deeper
    % eigenvalue indices are read off the rows themselves.
    if isempty(extras)
        timing_clause = [' Timing shows the time for the matrix assembly and ', ...
            'full spectrum computation using \texttt{MATLAB}''s ', ...
            '\lstinline{eig} function with the default settings.'];
    else
        timing_clause = '';
    end
    fprintf(fid, ['  \\caption{%s eigenvalues of Laplace operator on the %s ', ...
        'domain \\texttt{%s}. Numerical eigenvalues computed by each ', ...
        'discretisation %s with varying degrees of freedom (DOF)%s.%s}\n'], ...
        t.leading, cfg.pretty, latex_name(cfg.name), t.method_list, ...
        t.comparison_clause, timing_clause);
    % Distinct label: every layout may be \input into the same document.
    if isempty(extras)
        suffix = '';
    else
        suffix = '_extended';
    end
    fprintf(fid, '  \\label{tab:eigenvalues_head_%s_transposed%s}\n', cfg.name, suffix);
    fprintf(fid, '\\end{table}\n');
end


function write_eig_row(fid, rows, i, ndec)
%WRITE_EIG_ROW One eigenvalue row of the transposed table: index, then each run.
%
%   NDEC decimals per cell. A run holding fewer than i eigenvalues -- a grid
%   coarser than the index asked for, or a published MPS reference shorter than
%   the table -- gets "---", the marker the DOF and Time cells of the reference
%   column use.
    fprintf(fid, '    $\\lambda_{%d}$', i);
    for j = 1:numel(rows)
        row = rows(j);
        if i <= numel(row.eigs)
            s = sprintf('%.*f', ndec, row.eigs(i));
            if row.bold
                s = sprintf('\\textbf{%s}', s);
            end
        else
            s = '---';
        end
        fprintf(fid, ' & %s', s);
    end
    fprintf(fid, ' \\\\\n');
end


function s = index_list(idx)
%INDEX_LIST Eigenvalue indices as "$\lambda_{60}$, $\lambda_{80}$ and $\lambda_{100}$".
    parts = arrayfun(@(k) sprintf('$\\lambda_{%d}$', k), idx, 'UniformOutput', false);
    if numel(parts) == 1
        s = parts{1};
    else
        s = sprintf('%s and %s', strjoin(parts(1:end-1), ', '), parts{end});
    end
end


function s = latex_name(name)
%LATEX_NAME Escape underscores for \texttt{} in text mode.
    s = strrep(name, '_', '\_');
end
