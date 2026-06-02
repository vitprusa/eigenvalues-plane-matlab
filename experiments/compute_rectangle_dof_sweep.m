function compute_rectangle_dof_sweep()
%COMPUTE_RECTANGLE_DOF_SWEEP Generate the rectangle DOF-sweep spectra (if missing).
%
%   Computes the full Dirichlet-Laplacian spectrum on the rectangle
%   [0, 2*pi] x [0, pi] by DST, finite differences, finite elements (all dense
%   eig), and Chebyshev spectral collocation at several resolutions each, plus
%   the analytic ground truth lambda = m^2/4 + n^2. One CSV per (method,
%   setting) -- and one for the analytic spectrum -- is written into
%   results/eigenvalues_dof_sweep/. A CSV that already exists is skipped, so
%   re-running only fills in missing data.
%
%   Plot the result with PLOT_RECTANGLE_DOF_SWEEP.
%
%   Requires the PDE Toolbox (FEM) and Chebfun (Cheb).

    experiments_dir = fileparts(mfilename('fullpath'));
    project_root    = fileparts(experiments_dir);
    run(fullfile(project_root, 'startup.m'));

    out_dir = fullfile(project_root, 'results', 'eigenvalues_dof_sweep');
    if ~exist(out_dir, 'dir')
        mkdir(out_dir);
    end

    % Rectangle [0, 2*pi] x [0, pi] (fills its bounding box).
    box = [0 2*pi 0 pi];
    phi = @(x, y) indicator_rectangle(x, y, 0, 2*pi, 0, pi);   % DST / FD
    [x_range, y_range] = bounding_box(box(1), box(2), box(3), box(4));
    gd = [3; 4; 0; 2*pi; 2*pi; 0; 0; 0; pi; pi];               % FEM rectangle
    ns = char('R1')';
    sf = 'R1';

    % DOF settings (DST/FD include an extra ~3000-dof run at M = 77).
    M_grid   = [15 25 35 49 77];        % DST and FD (odd: M+1 even -> y=pi on grid)
    Hmax_fem = [0.35 0.25 0.18 0.13];   % FEM target mesh sizes
    N_cheb   = [10 20 30 40];           % Chebyshev orders

    % --- DST (full spectrum via dense eig) ---
    for k = 1:numel(M_grid)
        csv_file = fullfile(out_dir, sprintf('rectangle_dst_%d-eigenvalues.csv', k));
        if exist(csv_file, 'file'); continue; end
        M = M_grid(k);
        [L, info] = make_dst_laplace_mat_batched(x_range, y_range, M, phi);
        evals = sort(-real(eig(L)), 'ascend');
        write_sweep_csv(csv_file, evals, info.dofs, sprintf('M = %d', M), 'DST');
    end

    % --- FD (5-point stencil, dense eig) ---
    for k = 1:numel(M_grid)
        csv_file = fullfile(out_dir, sprintf('rectangle_fd_%d-eigenvalues.csv', k));
        if exist(csv_file, 'file'); continue; end
        M = M_grid(k);
        [evals, info] = fd_laplace_spectrum(struct('box', box, 'phi', phi, 'M', M));
        write_sweep_csv(csv_file, evals, info.dofs, sprintf('M = %d', M), 'FD');
    end

    % --- FEM (assembled K, M + dense eig) ---
    for k = 1:numel(Hmax_fem)
        csv_file = fullfile(out_dir, sprintf('rectangle_fem_%d-eigenvalues.csv', k));
        if exist(csv_file, 'file'); continue; end
        h = Hmax_fem(k);
        entry = struct('gd', gd, 'ns', ns, 'sf', sf, 'Hmax_eig', h, 'Hmax_solvepdeeig', h);
        [evals, info] = fem_laplace_spectrum(entry, "eig");
        write_sweep_csv(csv_file, evals, info.dofs, sprintf('Hmax = %g', h), 'FEM');
    end

    % --- Chebyshev spectral collocation ---
    for k = 1:numel(N_cheb)
        csv_file = fullfile(out_dir, sprintf('rectangle_cheb_%d-eigenvalues.csv', k));
        if exist(csv_file, 'file'); continue; end
        N = N_cheb(k);
        [evals, info] = chebfun_laplace_spectrum(box, N);
        write_sweep_csv(csv_file, evals, info.dofs, sprintf('N = %d', N), 'CHEB');
    end

    % --- Analytic ground truth: lambda = m^2/4 + n^2 (generous fixed length) ---
    csv_file = fullfile(out_dir, 'rectangle_analytic-eigenvalues.csv');
    if ~exist(csv_file, 'file')
        [Mm, Nn] = meshgrid(1:300, 1:150);
        v = sort(reshape((Mm.^2) / 4 + Nn.^2, [], 1), 'ascend');
        analytic = v(1:min(8000, numel(v)));
        fid = fopen(csv_file, 'w');
        fprintf(fid, '# Domain: rectangle (analytic ground truth, lambda = m^2/4 + n^2)\n');
        fprintf(fid, '# Computed %s\n', datestr(now, 'yyyy-mm-dd HH:MM:SS'));
        fclose(fid);
        n        = (1:numel(analytic))';
        lambda_n = analytic(:);
        writetable(table(n, lambda_n), csv_file, 'WriteMode', 'append', 'WriteVariableNames', true);
    end

    fprintf('Rectangle DOF-sweep data ready in %s\n', out_dir);
end


function write_sweep_csv(csv_file, evals, dofs, res, method)
    fid = fopen(csv_file, 'w');
    if fid == -1
        error('compute_rectangle_dof_sweep:cannotOpen', 'Could not open %s.', csv_file);
    end
    fprintf(fid, '# Domain: rectangle (%s, dense eig full spectrum)\n', method);
    fprintf(fid, '# Computed %s\n', datestr(now, 'yyyy-mm-dd HH:MM:SS'));
    fprintf(fid, '# Resolution %s, dofs = %d\n', res, dofs);
    fclose(fid);
    n        = (1:numel(evals))';
    lambda_n = evals(:);
    writetable(table(n, lambda_n), csv_file, 'WriteMode', 'append', 'WriteVariableNames', true);
end
