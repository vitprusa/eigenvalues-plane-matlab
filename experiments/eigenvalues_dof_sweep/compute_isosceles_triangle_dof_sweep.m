function compute_isosceles_triangle_dof_sweep()
%COMPUTE_ISOSCELES_TRIANGLE_DOF_SWEEP Generate the triangle DOF-sweep spectra (if missing).
%
%   Computes the full Dirichlet-Laplacian spectrum on the right isosceles
%   triangle with vertices (0,0), (pi,0), (pi,pi) (legs of length pi) by DST,
%   finite differences, and finite elements (all dense eig) at several
%   resolutions each (DST and FD include an extra ~3000-dof run at M = 77),
%   plus the analytic ground truth
%   lambda = m^2 + n^2 with integers m > n >= 1 (the antisymmetric modes of the
%   [0,pi]^2 square). One CSV per (method, setting) -- and one for the analytic
%   spectrum -- is written into results/eigenvalues_dof_sweep/. A CSV that
%   already exists is skipped, so re-running only fills in missing data.
%
%   Plot the result with PLOT_ISOSCELES_TRIANGLE_DOF_SWEEP.
%
%   Requires the PDE Toolbox (for the FEM part).

    dof_sweep_dir = fileparts(mfilename('fullpath'));
    project_root  = fileparts(fileparts(dof_sweep_dir));
    run(fullfile(project_root, 'startup.m'));

    out_dir = fullfile(project_root, 'results', 'eigenvalues_dof_sweep');
    if ~exist(out_dir, 'dir')
        mkdir(out_dir);
    end

    % Right isosceles triangle: 0 < x < pi, 0 < y < x (legs of length pi).
    box = [0 pi 0 pi];
    phi = @(x, y) indicator_isosceles_triangle(x, y, 0, pi, 0);   % DST / FD
    [x_range, y_range] = bounding_box(box(1), box(2), box(3), box(4));
    gd = [2; 3; 0; pi; pi; 0; 0; pi];     % FEM triangle (vertices (0,0),(pi,0),(pi,pi))
    ns = char('T1')';
    sf = 'T1';

    % DOF settings (DST/FD include an extra ~3000-dof run at M = 77).
    M_grid   = [15 25 35 50 77];        % DST and FD grid resolutions
    Hmax_fem = [0.20 0.13 0.09 0.06];   % FEM target mesh sizes

    % --- DST (full spectrum via dense eig) ---
    for k = 1:numel(M_grid)
        csv_file = fullfile(out_dir, sprintf('isosceles_triangle_dst_%d-eigenvalues.csv', k));
        if exist(csv_file, 'file'); continue; end
        M = M_grid(k);
        [L, info] = make_dst_laplace_mat_batched(x_range, y_range, M, phi);
        evals = sort(-real(eig(L)), 'ascend');
        write_sweep_csv(csv_file, evals, info.dofs, sprintf('M = %d', M), 'DST');
    end

    % --- FD (5-point stencil, dense eig) ---
    for k = 1:numel(M_grid)
        csv_file = fullfile(out_dir, sprintf('isosceles_triangle_fd_%d-eigenvalues.csv', k));
        if exist(csv_file, 'file'); continue; end
        M = M_grid(k);
        [evals, info] = fd_laplace_spectrum(struct('box', box, 'phi', phi, 'M', M));
        write_sweep_csv(csv_file, evals, info.dofs, sprintf('M = %d', M), 'FD');
    end

    % --- FEM (assembled K, M + dense eig) ---
    for k = 1:numel(Hmax_fem)
        csv_file = fullfile(out_dir, sprintf('isosceles_triangle_fem_%d-eigenvalues.csv', k));
        if exist(csv_file, 'file'); continue; end
        h = Hmax_fem(k);
        entry = struct('gd', gd, 'ns', ns, 'sf', sf, 'Hmax_eig', h, 'Hmax_solvepdeeig', h);
        [evals, info] = fem_laplace_spectrum(entry, "eig");
        write_sweep_csv(csv_file, evals, info.dofs, sprintf('Hmax = %g', h), 'FEM');
    end

    % --- Analytic ground truth: lambda = m^2 + n^2, m > n >= 1 (legs pi) ---
    csv_file = fullfile(out_dir, 'isosceles_triangle_analytic-eigenvalues.csv');
    if ~exist(csv_file, 'file')
        [Mm, Nn] = meshgrid(1:200, 1:200);
        keep = Mm > Nn;                                % m > n >= 1
        v = sort(Mm(keep).^2 + Nn(keep).^2, 'ascend');
        analytic = v(1:min(8000, numel(v)));
        fid = fopen(csv_file, 'w');
        fprintf(fid, '# Domain: isosceles_triangle (analytic ground truth, lambda = m^2 + n^2, m > n >= 1)\n');
        fprintf(fid, '# Computed %s\n', datestr(now, 'yyyy-mm-dd HH:MM:SS'));
        fclose(fid);
        n        = (1:numel(analytic))';
        lambda_n = analytic(:);
        writetable(table(n, lambda_n), csv_file, 'WriteMode', 'append', 'WriteVariableNames', true);
    end

    fprintf('Isosceles-triangle DOF-sweep data ready in %s\n', out_dir);
end


function write_sweep_csv(csv_file, evals, dofs, res, method)
    fid = fopen(csv_file, 'w');
    if fid == -1
        error('compute_isosceles_triangle_dof_sweep:cannotOpen', 'Could not open %s.', csv_file);
    end
    fprintf(fid, '# Domain: isosceles_triangle (%s, dense eig full spectrum)\n', method);
    fprintf(fid, '# Computed %s\n', datestr(now, 'yyyy-mm-dd HH:MM:SS'));
    fprintf(fid, '# Resolution %s, dofs = %d\n', res, dofs);
    fclose(fid);
    n        = (1:numel(evals))';
    lambda_n = evals(:);
    writetable(table(n, lambda_n), csv_file, 'WriteMode', 'append', 'WriteVariableNames', true);
end
