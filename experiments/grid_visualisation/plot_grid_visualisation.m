function plot_grid_visualisation(name)
%PLOT_GRID_VISUALISATION Visualise the DST grid and domain mask for each domain.
%
%   plot_grid_visualisation() draws one figure per domain in DOMAIN_CATALOG_DST,
%   using the SAME bounding box, grid spacing, and mask as the DST "full"
%   spectrum (resolution M_full), and writes a EPS into
%   results/grid_visualisation/.
%
%   plot_grid_visualisation(name) restricts the run to the single domain "name"
%   (e.g. "L_shaped"). Pass "" or [] to keep all domains.
%
%   Each figure shows:
%     - the continuous domain, shaded grey (the region where phi > 0);
%     - every grid point of the bounding-box grid as an open black circle;
%     - the masked grid points -- the interior degrees of freedom that hold the
%       function values in the DST discretisation -- as filled black circles.
%
%   The grid and mask are taken from MAKE_DST_LAPLACE_OP so that the picture is
%   exactly the grid the DST method solves on.
%
%   See also DOMAIN_CATALOG_DST, MAKE_DST_LAPLACE_OP, COMPUTE_SPECTRUM_DST.

    if nargin < 1
        name = '';
    end

    gv_dir          = fileparts(mfilename('fullpath'));
    experiments_dir = fileparts(gv_dir);
    project_root    = fileparts(experiments_dir);
    run(fullfile(project_root, 'startup.m'));
    addpath(experiments_dir);   % domain_catalog_dst

    out_dir = fullfile(project_root, 'results', 'grid_visualisation');
    if ~exist(out_dir, 'dir')
        mkdir(out_dir);
    end

    cases = domain_catalog_dst(name);
    for i = 1:numel(cases)
        fprintf('=== %s : grid visualisation (M_full = %d) ===\n', ...
            cases(i).name, cases(i).M_full);
        plot_one(cases(i), out_dir);
    end
end


function plot_one(c, out_dir)
    [x_range, y_range] = bounding_box(c.box(1), c.box(2), c.box(3), c.box(4));
    M = c.M_full;

    % Reuse the DST factory for the exact grid and mask it solves on.
    [~, info] = make_dst_laplace_op(x_range, y_range, M, c.phi);
    X    = info.X;
    Y    = info.Y;
    mask = info.global_mask;

    % Render text (title, legend) in the standard LaTeX Computer Modern font.
    fig = figure('Visible', 'off', 'Position', [100 100 760 760], ...
        'defaultAxesTickLabelInterpreter', 'latex', ...
        'defaultTextInterpreter',          'latex', ...
        'defaultLegendInterpreter',        'latex');
    ax = axes(fig);
    hold(ax, 'on');

    % --- Domain, shaded grey ---------------------------------------------
    % Sample the indicator on a fine grid and fill the region phi > 0 grey.
    nf = 800;
    xf = linspace(x_range(1), x_range(2), nf);
    yf = linspace(y_range(1), y_range(2), nf);
    [Xf, Yf] = meshgrid(xf, yf);
    dom = c.phi(Xf, Yf) > 0;
    inside = double(dom);
    inside(~dom) = NaN;                 % leave the exterior unpainted
    pc = pcolor(ax, Xf, Yf, inside);
    set(pc, 'EdgeColor', 'none', 'HandleVisibility', 'off');
    colormap(ax, [0.82 0.82 0.82]);     % single grey
    clim(ax, [0 1]);

    % Domain boundary as a thin line (the 0.5 level set of the indicator).
    contour(ax, Xf, Yf, double(dom), [0.5 0.5], 'LineColor', 'k', ...
        'LineWidth', 0.75, 'HandleVisibility', 'off');

    % Proxy handles so the grey domain and its boundary get clean legend entries.
    hdom = patch(ax, NaN, NaN, [0.82 0.82 0.82], 'EdgeColor', 'none', ...
        'DisplayName', 'domain ($\varphi > 0$)');
    hbnd = plot(ax, NaN, NaN, '-', 'Color', 'k', 'LineWidth', 0.75, ...
        'DisplayName', 'domain boundary');

    % --- Grid points ------------------------------------------------------
    hgrid = plot(ax, X(:), Y(:), 'o', 'MarkerEdgeColor', 'k', ...
        'MarkerFaceColor', 'none', 'MarkerSize', 4, 'LineStyle', 'none', ...
        'DisplayName', 'grid point');
    hmask = plot(ax, X(mask), Y(mask), 'o', 'MarkerEdgeColor', 'k', ...
        'MarkerFaceColor', 'k', 'MarkerSize', 4, 'LineStyle', 'none', ...
        'DisplayName', 'masked grid point ($\mathtt{DOF}$)');

    hold(ax, 'off');
    axis(ax, 'equal');
    margin = info.h;
    xlim(ax, [x_range(1) - margin, x_range(2) + margin]);
    ylim(ax, [y_range(1) - margin, y_range(2) + margin]);
    box(ax, 'on');

    pretty = strrep(c.name, '_', '\_');
    title(ax, sprintf('%s --- DST grid ($M = %d$, dofs $= %d$)', ...
        pretty, M, info.dofs));
    legend(ax, [hdom hbnd hgrid hmask], 'Location', 'northeastoutside', ...
        'FontSize', 9, 'Box', 'off');

    eps_file = fullfile(out_dir, sprintf('%s_grid.eps', c.name));
    try
        exportgraphics(fig, eps_file, 'ContentType', 'vector');
    catch
        print(fig, eps_file, '-depsc2', '-painters');
    end
    close(fig);
    fprintf('Wrote %s\n', eps_file);
end
