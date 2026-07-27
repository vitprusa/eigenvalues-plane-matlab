function eigenvalues_rectangle_tables_paper(M, margin)
%EIGENVALUES_RECTANGLE_TABLES_PAPER LaTeX snippet with the rectangle match and order tables.
%
%   eigenvalues_rectangle_tables_paper()        uses M = 9 and a default margin.
%   eigenvalues_rectangle_tables_paper(M)       uses resolution M.
%   eigenvalues_rectangle_tables_paper(M, margin) extends the analytic (m, n)
%                                         table margin indices past the resolved
%                                         modes in each direction.
%
%   Paper counterpart of EIGENVALUES_RECTANGLE_MATCH and
%   EIGENVALUES_RECTANGLE_ORDER combined: it performs the SAME computation as
%   those two scripts -- the full DST-Laplacian spectrum of the 2:1 catalog
%   rectangle [0, 2*pi] x [0, pi], the analytic lattice
%
%       lambda_{m,n} = (m*pi/Lx)^2 + (n*pi/Ly)^2 = m^2/4 + n^2,   m, n = 1, 2, ...
%
%   and the greedy multiplicity-aware match of the computed spectrum to the
%   resolved Mx-by-My block -- but instead of two markdown files it writes a
%   single LaTeX snippet holding BOTH tables:
%
%     1. values   -- the analytic eigenvalues lambda_{m,n}, and
%     2. ordering -- the position (rank) of each lambda_{m,n} in the
%                    magnitude-ordered sequence of all analytic eigenvalues.
%
%   In both tables the cells the grid resolves (matched to a computed eigenvalue
%   to round-off) are set in bold and the unresolved border in normal weight.
%   The snippet goes to results_paper/eigenvalues_indexing/rectangle_tables_M<M>.tex.
%
%   The default M = 9 reproduces the worked example in the markdown results
%   (Mx = 9, My = 4, table m = 1..14, n = 1..9).
%
%   See also EIGENVALUES_RECTANGLE_MATCH, EIGENVALUES_RECTANGLE_ORDER.

    if nargin < 1 || isempty(M)
        M = 9;
    end
    if nargin < 2 || isempty(margin)
        margin = 5;     % table reaches Mx+margin = 14 rows, My+margin = 9 cols
    end

    % This file lives in experiments_paper/eigenvalues_indexing_paper/; put the
    % project sources on the path.
    script_dir      = fileparts(mfilename('fullpath'));
    experiments_dir = fileparts(script_dir);
    project_root    = fileparts(experiments_dir);
    run(fullfile(project_root, 'startup.m'));

    % 2:1 catalog rectangle [0, Lx] x [0, Ly]. Analytic spectrum is
    % lambda_{m,n} = (m*pi/Lx)^2 + (n*pi/Ly)^2 = m^2/4 + n^2.
    Lx = 2 * pi;
    Ly = pi;

    % --- 1. Computed full spectrum at resolution M -----------------------
    [x_range, y_range] = bounding_box(0, Lx, 0, Ly);
    phi = @(x, y) indicator_rectangle(x, y, 0, Lx, 0, Ly);
    [L, info] = make_dst_laplace_mat_batched(x_range, y_range, M, phi);
    lambda_num = sort(-real(eig(L)), 'ascend');   % positive, ascending

    % Modes resolved per direction = interior grid points per direction.
    Mx = numel(info.x_vec) - 2;
    My = numel(info.y_vec) - 2;
    if Mx * My ~= info.dofs
        error('eigenvalues_rectangle_tables_paper:misaligned', ...
            ['Grid does not cleanly resolve the rectangle (Mx*My = %d, ', ...
             'dofs = %d). Pick M so both edges land on grid lines ', ...
             '(M odd for the 2:1 rectangle).'], Mx * My, info.dofs);
    end

    fprintf('Rectangle [0, %.4f] x [0, %.4f], M = %d, dofs = %d (Mx = %d, My = %d)\n', ...
        Lx, Ly, M, info.dofs, Mx, My);

    % --- 2. Analytic eigenvalue lattice indexed by (m, n) ---------------
    mm = (1:(Mx + margin))';
    nn = (1:(My + margin));
    lambda_tab = (mm * pi / Lx).^2 + (nn * pi / Ly).^2;

    % --- 3. Match computed eigenvalues to lattice cells -----------------
    % Greedy, multiplicity-aware matching by value; resolved cells (m <= Mx,
    % n <= My) are visited first so they claim their computed eigenvalue ahead
    % of any coincident out-of-grid cell of equal value.
    cell_val = lambda_tab(:);
    in_grid  = false(numel(mm), numel(nn));
    in_grid(1:Mx, 1:My) = true;
    order    = [find(in_grid(:)); find(~in_grid(:))];

    reltol  = 1e-6;
    used    = false(numel(lambda_num), 1);
    matched = false(numel(mm), numel(nn));
    for c = order'
        target = cell_val(c);
        cand = find(~used & abs(lambda_num - target) <= reltol * target, 1);
        if ~isempty(cand)
            used(cand)  = true;
            matched(c)  = true;
        end
    end
    n_matched   = nnz(matched);
    n_unmatched = nnz(~used);

    % --- 4. Position (rank) of each eigenvalue in the ordered list ------
    rank_tab = eigenvalue_positions(lambda_tab, Lx, Ly);

    if n_unmatched > 0
        warning('eigenvalues_rectangle_tables_paper:unmatched', ...
            ['%d computed eigenvalues did not match any (m, n) cell; ', ...
             'increase margin so the lattice covers all resolved modes.'], ...
            n_unmatched);
    end

    meta = struct('Lx', Lx, 'Ly', Ly, 'M', M, 'dofs', info.dofs, ...
        'Mx', Mx, 'My', My, 'n_matched', n_matched, ...
        'n_cells', numel(matched), 'n_computed', numel(lambda_num), ...
        'n_unmatched', n_unmatched);

    % --- 5. Write the LaTeX snippet with both tables --------------------
    out_dir = fullfile(project_root, 'results_paper', 'eigenvalues_indexing');
    if ~exist(out_dir, 'dir')
        mkdir(out_dir);
    end
    tex_file = fullfile(out_dir, sprintf('rectangle_tables_M%d.tex', M));
    write_latex(tex_file, lambda_tab, rank_tab, matched, mm, nn, meta);
    fprintf('Wrote %s\n', tex_file);
end


function rank_tab = eigenvalue_positions(lambda_tab, Lx, Ly)
%EIGENVALUE_POSITIONS Position of each lattice cell in the magnitude-ordered spectrum.
%   Returns a matrix the size of lambda_tab whose (r, c) entry is the rank of
%   lambda_{m,n} (1 = smallest) among ALL analytic eigenvalues, not just those
%   in the table. The ranking lattice is grown until it contains every (m, n)
%   whose eigenvalue does not exceed the largest one displayed.
    [nrow, ncol] = size(lambda_tab);
    thr = lambda_tab(nrow, ncol);                 % largest displayed eigenvalue
    Mrank = max(nrow, ceil((Lx / pi) * sqrt(thr)) + 1);
    Nrank = max(ncol, ceil((Ly / pi) * sqrt(thr)) + 1);

    [Mr, Nr] = ndgrid(1:Mrank, 1:Nrank);
    val = (Mr * pi / Lx).^2 + (Nr * pi / Ly).^2;

    % Sort by magnitude; round first so genuinely degenerate eigenvalues tie
    % exactly (the smallest spacing here is 1/4), then break ties by (m, n) for
    % a deterministic sequence. Position = index in that sequence.
    keys = [round(val(:), 9), Mr(:), Nr(:)];
    [~, ord] = sortrows(keys);
    pos = zeros(Mrank * Nrank, 1);
    pos(ord) = 1:numel(ord);
    pos = reshape(pos, Mrank, Nrank);

    rank_tab = pos(1:nrow, 1:ncol);
end


function write_latex(tex_file, lambda_tab, rank_tab, matched, mm, nn, meta)
%WRITE_LATEX Write both tables (values and ordering) as one subfig LaTeX snippet.
    fid = fopen(tex_file, 'w');
    if fid == -1
        error('eigenvalues_rectangle_tables_paper:cannotOpen', ...
            'Could not open %s for writing.', tex_file);
    end
    closer = onCleanup(@() fclose(fid));

    % --- Header comment block, matching the house style of the other snippets.
    fprintf(fid, ['%% Rectangle DST eigenvalue indexing tables (values and ', ...
        'ordering), M = %d.\n'], meta.M);
    fprintf(fid, '%%\n');
    fprintf(fid, '%% Requires in the preamble:\n');
    fprintf(fid, '%%   \\usepackage{booktabs}\n');
    fprintf(fid, '%%   \\usepackage{subfig}\n');
    fprintf(fid, '%%\n');
    fprintf(fid, ['%% The full DST-Laplacian spectrum of the 2:1 rectangle ', ...
        '[0, 2pi] x [0, pi] at\n']);
    fprintf(fid, ['%% resolution M = %d (dofs = %d) is matched against the ', ...
        'analytic Dirichlet\n'], meta.M, meta.dofs);
    fprintf(fid, ['%% eigenvalues lambda_{m,n} = m^2/4 + n^2. The grid ', ...
        'resolves Mx = %d modes in x\n'], meta.Mx);
    fprintf(fid, ['%% and My = %d modes in y; matched (resolved) cells are ', ...
        'set in bold.\n'], meta.My);
    fprintf(fid, '\n');

    % booktabs horizontal rules (\toprule etc.) plus one vertical rule after the
    % first column, spanning the whole table (header and body).
    colspec = ['r|', repmat('r', 1, numel(nn))];

    % --- The float with two stacked sub-tables.
    fprintf(fid, '\\begin{table}\n');
    fprintf(fid, '  \\centering\n');

    % Sub-table (a): the analytic eigenvalues.
    fprintf(fid, ['  \\subfloat[Computed eigenvalues in \\textbf{bold}. The ', ...
        'computational grid resolves $M_x = %d$ modes in the horizontal ', ...
        'direction $x$ and $M_y = %d$ modes in the vertical direction $y$, ', ...
        'total number of computed eigenvalues is $M_x \\times M_y = ', ...
        '\\text{DOF}$. Computed eigenvalues match---up to round-off error---', ...
        'the corresponding analytic eigenvalues, which are computed for all ', ...
        'pairs $(m, n)$ and are shown in normal weight.]{%%\n'], ...
        meta.Mx, meta.My);
    write_tabular(fid, lambda_tab, matched, mm, nn, colspec, @(v) sprintf('%.2f', v));
    fprintf(fid, '\n  }\n');
    fprintf(fid, '  \\\\\n');

    % Sub-table (b): the positions.
    fprintf(fid, ['  \\subfloat[Sequential position of analytical eigenvalues ', ...
        '$\\lambda_{m,n}$ in the magnitude-ordered sequence of all analytic ', ...
        'eigenvalues.]{%%\n']);
    write_tabular(fid, rank_tab, matched, mm, nn, colspec, @(v) sprintf('%d', v));
    fprintf(fid, '\n  }\n');

    fprintf(fid, ['  \\caption{Analytical/numerical eigenvalues of Laplace ', ...
        'operator on the rectangular domain $[0, %s] \\times [0, %s]$, zero ', ...
        'Dirichlet boundary condition. The analytical eigenvalues of the ', ...
        'continuous eigenvalue problem are given by the formula~', ...
        '$\\lambda_{m,n} = \\frac{m^2}{4} + n^2$. (Eigenvalues indexed by~', ...
        '$(m, n)$, rows are indexed by~$m$, columns by~$n$.) Numerical ', ...
        'eigenvalues computed using the \\texttt{DST} based discretisation, ', ...
        '$M = %d$, $\\text{DOF} = %d$.}\n'], ...
        pi_label(meta.Lx), pi_label(meta.Ly), meta.M, meta.dofs);
    fprintf(fid, '  \\label{tab:rectangle_indexing_M%d}\n', meta.M);
    fprintf(fid, '\\end{table}\n');
end


function write_tabular(fid, data, matched, mm, nn, colspec, fmt)
%WRITE_TABULAR Emit one tabular for DATA indexed by m (rows) and n (cols).
%   Matched cells are wrapped in \textbf{...}. fmt maps a value to its string.
    fprintf(fid, '    \\begin{tabular}{%s}\n', colspec);
    fprintf(fid, '      \\toprule\n');

    % Header row: n indices across the top.
    fprintf(fid, '      $m \\backslash n$');
    for n = nn
        fprintf(fid, ' & %d', n);
    end
    fprintf(fid, ' \\\\\n');
    fprintf(fid, '      \\midrule\n');

    % One row per m.
    for r = 1:numel(mm)
        fprintf(fid, '      %d', mm(r));
        for c = 1:numel(nn)
            s = fmt(data(r, c));
            if matched(r, c)
                fprintf(fid, ' & \\textbf{%s}', s);
            else
                fprintf(fid, ' & %s', s);
            end
        end
        fprintf(fid, ' \\\\\n');
    end

    fprintf(fid, '      \\bottomrule\n');
    fprintf(fid, '    \\end{tabular}');
end


function s = pi_label(v)
%PI_LABEL LaTeX label for a length that is an integer multiple of pi.
    k = v / pi;
    if abs(k - round(k)) < 1e-9
        k = round(k);
        if k == 1
            s = '\pi';
        else
            s = sprintf('%d\\pi', k);
        end
    else
        s = sprintf('%.4f', v);
    end
end
