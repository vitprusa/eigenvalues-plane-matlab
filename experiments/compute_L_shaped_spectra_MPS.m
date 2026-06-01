function compute_L_shaped_spectra_MPS()
%COMPUTE_L_SHAPED_SPECTRA_MPS Run the MPS L-shaped eigenvalue experiment.
%
%   compute_L_shaped_spectra_MPS() computes the leading Dirichlet-Laplacian
%   eigenvalues on the L-shaped region by the Method of Particular Solutions
%   and writes them to results/mps/L_shaped_eigenvalues_MPS.csv with a
%   metadata header (columns n, lambda_n).
%
%   The computation itself is the Betcke & Trefethen "Ldrum" code: this runner
%   executes the core script src/mps/Ldrum_modified.m (which leaves the
%   eigenvalues in evals and the parameters K, N, np, lammax in the
%   workspace), then exports the results. To change the number of eigenvalues
%   or the accuracy, edit src/mps/Ldrum_modified.m.
%
%   See also LDRUM_MODIFIED, LDRUM.

    % This file lives in experiments/; put the project sources on the path.
    project_root = fileparts(fileparts(mfilename('fullpath')));
    run(fullfile(project_root, 'startup.m'));

    % Run the core MPS script (on the path via src/mps). It executes in this
    % workspace, leaving evals and the parameters K, N, np, lammax available.
    Ldrum_modified;

    out_dir = fullfile(project_root, 'results', 'mps');
    if ~exist(out_dir, 'dir')
        mkdir(out_dir);
    end
    csv_file = fullfile(out_dir, 'L_shaped_eigenvalues_MPS.csv');

    fid = fopen(csv_file, 'w');
    if fid == -1
        error('compute_L_shaped_spectra_MPS:cannotOpen', 'Could not open %s for writing.', csv_file);
    end
    fprintf(fid, '# Domain: L_shaped (MPS, method of particular solutions)\n');
    fprintf(fid, '# Computed %s\n', datestr(now, 'yyyy-mm-dd HH:MM:SS'));
    fprintf(fid, '# Betcke & Trefethen, SIAM Review 47(3):469-491, 2005\n');
    fprintf(fid, '# K = %d, accuracy N = %d, points np = %d, lammax = %g\n', K, N, np, lammax);
    fclose(fid);

    n        = (1:numel(evals))';
    lambda_n = evals(:);
    writetable(table(n, lambda_n), csv_file, 'WriteMode', 'append', 'WriteVariableNames', true);

    fprintf('Wrote %d eigenvalues to %s\n', numel(evals), csv_file);
end
