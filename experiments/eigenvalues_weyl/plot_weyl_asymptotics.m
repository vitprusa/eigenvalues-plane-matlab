function plot_weyl_asymptotics(name)
%PLOT_WEYL_ASYMPTOTICS Visual Weyl-asymptotics check for each DST-catalog domain.
%
%   plot_weyl_asymptotics() processes every domain in DOMAIN_CATALOG_DST: it
%   computes the full Dirichlet-Laplacian spectrum at the "full" resolution
%   M_full (dense eig on MAKE_DST_LAPLACE_MAT_BATCHED), computes the domain
%   area, and draws the Weyl check
%
%       lambda_n / n  -->  4*pi / Area      as n grows,
%
%   the two-dimensional Weyl law N(lambda) ~ (Area/(4*pi)) lambda written per
%   eigenvalue. One EPS per domain is written into results/eigenvalues_weyl/,
%   together with a summary table weyl_areas.md.
%
%   plot_weyl_asymptotics(name) restricts the run to the single domain "name"
%   (e.g. "L_shaped"). Pass "" or [] to keep all domains.
%
%   The area is computed by fine-grid quadrature of the indicator (the
%   fraction of the bounding box with phi > 0, times the box area), so it is
%   the continuous geometric area, independent of the DST resolution.
%
%   This is the refactored, catalog-driven counterpart of the asymptotic
%   blocks in the old dst_laplace_grid_eigenvalues_*_test scripts.
%
%   See also DOMAIN_CATALOG_DST, MAKE_DST_LAPLACE_MAT_BATCHED.

    if nargin < 1
        name = '';
    end

    weyl_dir        = fileparts(mfilename('fullpath'));
    experiments_dir = fileparts(weyl_dir);
    project_root    = fileparts(experiments_dir);
    run(fullfile(project_root, 'startup.m'));
    addpath(experiments_dir);   % domain_catalog_dst

    out_dir = fullfile(project_root, 'results', 'eigenvalues_weyl');
    if ~exist(out_dir, 'dir')
        mkdir(out_dir);
    end

    cases = domain_catalog_dst(name);
    rows  = struct('name', {}, 'area', {}, 'dofs', {}, 'C', {});
    for i = 1:numel(cases)
        fprintf('=== %s : Weyl asymptotics (M_full = %d) ===\n', ...
            cases(i).name, cases(i).M_full);
        rows(end+1) = plot_one(cases(i), out_dir); %#ok<AGROW>
    end

    write_summary(fullfile(out_dir, 'weyl_areas.md'), rows);
end


function row = plot_one(c, out_dir)
    [x_range, y_range] = bounding_box(c.box(1), c.box(2), c.box(3), c.box(4));
    M = c.M_full;

    % Full Dirichlet-Laplacian spectrum at the "full" resolution.
    [L, info] = make_dst_laplace_mat_batched(x_range, y_range, M, c.phi);
    % lambda = sort(-real(eig(L)), 'ascend');     % positive eigenvalues, ascending
    % positive eigenvalues, ascending
    lambda = sort(-eig(dst_laplace_symmetrise(L)), 'ascend');
    dofs   = info.dofs;
    n      = (1:dofs)';

    % Continuous geometric area by fine-grid quadrature of the indicator.
    area = indicator_area(c.phi, x_range, y_range);

    % Weyl law per eigenvalue: lambda_n / n -> 4*pi / Area.
    E = lambda ./ n;
    C = 4 * pi / area;

    % Render text in the standard LaTeX Computer Modern font.
    fig = figure('Visible', 'off', 'Position', [100 100 900 620], ...
        'defaultAxesTickLabelInterpreter', 'latex', ...
        'defaultTextInterpreter',          'latex', ...
        'defaultLegendInterpreter',        'latex');
    ax = axes(fig);
    hold(ax, 'on');

    plot(ax, n, E, 'o', 'MarkerEdgeColor', [0 0.45 0.74], 'MarkerSize', 3, ...
        'LineStyle', 'none', 'DisplayName', '$\lambda_n / n$');
    yline(ax, C, 'k--', 'LineWidth', 1.2, ...
        'DisplayName', sprintf('Weyl limit $4\\pi/A = %.4f$', C));

    hold(ax, 'off');
    grid(ax, 'on'); box(ax, 'on');
    xlim(ax, [0 dofs]);
    xlabel(ax, 'eigenvalue index $n$');
    ylabel(ax, '$\lambda_n / n$');
    pretty = strrep(c.name, '_', '\_');
    title(ax, sprintf('%s --- Weyl asymptotics ($M = %d$, dofs $= %d$, $A = %.4f$)', ...
        pretty, M, dofs, area));
    legend(ax, 'Location', 'northwest', 'FontSize', 9, 'Box', 'off');

    eps_file = fullfile(out_dir, sprintf('%s_weyl.eps', c.name));
    try
        exportgraphics(fig, eps_file, 'ContentType', 'vector');
    catch
        print(fig, eps_file, '-depsc2', '-painters');
    end
    close(fig);
    fprintf('Wrote %s  (area = %.6f, 4*pi/A = %.6f)\n', eps_file, area, C);

    row = struct('name', c.name, 'area', area, 'dofs', dofs, 'C', C);
end


function area = indicator_area(phi, x_range, y_range)
%INDICATOR_AREA Continuous area of {phi > 0} by fine-grid quadrature.
    nf = 4000;
    xf = linspace(x_range(1), x_range(2), nf);
    yf = linspace(y_range(1), y_range(2), nf);
    [Xf, Yf] = meshgrid(xf, yf);
    inside = phi(Xf, Yf) > 0;
    box_area = (x_range(2) - x_range(1)) * (y_range(2) - y_range(1));
    area = mean(inside(:)) * box_area;
end


function write_summary(md_file, rows)
%WRITE_SUMMARY Write a markdown table of the per-domain Weyl quantities.
    fid = fopen(md_file, 'w');
    if fid == -1
        error('plot_weyl_asymptotics:cannotOpen', 'Could not open %s.', md_file);
    end
    fprintf(fid, '# Weyl asymptotics --- domain areas and limits\n\n');
    fprintf(fid, ['Per-domain area (fine-grid quadrature of the indicator), ', ...
                  'degrees of freedom at the full resolution, and the Weyl ', ...
                  'limit $4\\pi/A$ that $\\lambda_n/n$ approaches.\n\n']);
    fprintf(fid, '| Domain | Area $A$ | dofs | $4\\pi/A$ |\n');
    fprintf(fid, '|--------|---------:|-----:|---------:|\n');
    for i = 1:numel(rows)
        fprintf(fid, '| %s | %.4f | %d | %.4f |\n', ...
            strrep(rows(i).name, '_', '\_'), rows(i).area, rows(i).dofs, rows(i).C);
    end
    fclose(fid);
    fprintf('Wrote %s\n', md_file);
end
