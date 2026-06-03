function compute_L_shaped_dof_sweep()
%COMPUTE_L_SHAPED_DOF_SWEEP Generate the L-shaped DOF-sweep spectra (if missing).
%
%   Computes the full Dirichlet-Laplacian spectrum on the L-shaped domain by
%   DST, finite differences, and finite elements (all dense eig) at four
%   resolutions (degree-of-freedom settings) each, and writes one CSV per
%   (method, setting) into results/eigenvalues_dof_sweep/. A setting whose CSV
%   already exists is skipped, so re-running only fills in missing data.
%
%   Plot the result with PLOT_L_SHAPED_DOF_SWEEP.
%
%   Requires the PDE Toolbox (for the FEM part).

    dof_sweep_dir = fileparts(mfilename('fullpath'));
    project_root  = fileparts(fileparts(dof_sweep_dir));
    run(fullfile(project_root, 'startup.m'));

    out_dir = fullfile(project_root, 'results', 'eigenvalues_dof_sweep');
    if ~exist(out_dir, 'dir')
        mkdir(out_dir);
    end

    % L-shaped domain: [-1,1]^2 minus the upper-right unit square [0,1]^2.
    box = [-1 1 -1 1];
    phi = @(x, y) indicator_L_shaped(x, y, -1, 0, 1, -1, 0, 1);   % DST / FD
    [x_range, y_range] = bounding_box(box(1), box(2), box(3), box(4));
    gd = [[3; 4; -1; 1; 1; -1; -1; -1; 1; 1], [3; 4; 0; 1; 1; 0; 0; 0; 1; 1]];   % FEM
    ns = char('R1', 'R2')';
    sf = 'R1-R2';

    M_grid   = [15 25 35 50];           % DST and FD grid resolutions
    Hmax_fem = [0.20 0.13 0.09 0.06];   % FEM target mesh sizes

    % --- DST (full spectrum via dense eig) ---
    for k = 1:numel(M_grid)
        csv_file = fullfile(out_dir, sprintf('L_shaped_dst_%d-eigenvalues.csv', k));
        if exist(csv_file, 'file'); continue; end
        M = M_grid(k);
        [L, info] = make_dst_laplace_mat_batched(x_range, y_range, M, phi);
        evals = sort(-real(eig(L)), 'ascend');
        write_sweep_csv(csv_file, evals, info.dofs, sprintf('M = %d', M), 'DST');
    end

    % --- FD (5-point stencil, dense eig) ---
    for k = 1:numel(M_grid)
        csv_file = fullfile(out_dir, sprintf('L_shaped_fd_%d-eigenvalues.csv', k));
        if exist(csv_file, 'file'); continue; end
        M = M_grid(k);
        [evals, info] = fd_laplace_spectrum(struct('box', box, 'phi', phi, 'M', M));
        write_sweep_csv(csv_file, evals, info.dofs, sprintf('M = %d', M), 'FD');
    end

    % --- FEM (assembled K, M + dense eig) ---
    for k = 1:numel(Hmax_fem)
        csv_file = fullfile(out_dir, sprintf('L_shaped_fem_%d-eigenvalues.csv', k));
        if exist(csv_file, 'file'); continue; end
        h = Hmax_fem(k);
        entry = struct('gd', gd, 'ns', ns, 'sf', sf, 'Hmax_eig', h, 'Hmax_solvepdeeig', h);
        [evals, info] = fem_laplace_spectrum(entry, "eig");
        write_sweep_csv(csv_file, evals, info.dofs, sprintf('Hmax = %g', h), 'FEM');
    end

    fprintf('L-shaped DOF-sweep data ready in %s\n', out_dir);
end


function write_sweep_csv(csv_file, evals, dofs, res, method)
    fid = fopen(csv_file, 'w');
    if fid == -1
        error('compute_L_shaped_dof_sweep:cannotOpen', 'Could not open %s.', csv_file);
    end
    fprintf(fid, '# Domain: L_shaped (%s, dense eig full spectrum)\n', method);
    fprintf(fid, '# Computed %s\n', datestr(now, 'yyyy-mm-dd HH:MM:SS'));
    fprintf(fid, '# Resolution %s, dofs = %d\n', res, dofs);
    fclose(fid);
    n        = (1:numel(evals))';
    lambda_n = evals(:);
    writetable(table(n, lambda_n), csv_file, 'WriteMode', 'append', 'WriteVariableNames', true);
end
