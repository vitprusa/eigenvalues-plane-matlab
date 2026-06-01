function compute_rectangle_spectra_cheb()
%COMPUTE_RECTANGLE_SPECTRA_CHEB Run the Chebfun rectangle eigenvalue experiment.
%
%   compute_rectangle_spectra_cheb() computes the Dirichlet-Laplacian spectrum
%   on the rectangle [0,pi]x[0,pi] by Chebyshev spectral collocation (Chebfun
%   diffmat) and writes it to results/cheb/rectangle_eigenvalues_cheb.csv with
%   a metadata header (columns n, lambda_n).
%
%   The computation is the core script
%   src/cheb/chebfun_eigenvalues_rectangle.m, which this runner executes
%   (leaving the eigenvalues in eigs_chebfun and the parameters a, b, c, d, N
%   in the workspace); this file owns the path setup and the CSV export. To
%   change the resolution N or the box, edit that script. Requires Chebfun.
%
%   See also CHEBFUN_EIGENVALUES_RECTANGLE.

    % This file lives in experiments/; put the project sources on the path.
    project_root = fileparts(fileparts(mfilename('fullpath')));
    run(fullfile(project_root, 'startup.m'));

    % Run the core script (on the path via src/cheb). It executes in this
    % workspace, leaving eigs_chebfun and the parameters a, b, c, d, N.
    chebfun_eigenvalues_rectangle;

    out_dir = fullfile(project_root, 'results', 'cheb');
    if ~exist(out_dir, 'dir')
        mkdir(out_dir);
    end
    csv_file = fullfile(out_dir, 'rectangle_eigenvalues_cheb.csv');

    fid = fopen(csv_file, 'w');
    if fid == -1
        error('compute_rectangle_spectra_cheb:cannotOpen', 'Could not open %s for writing.', csv_file);
    end
    fprintf(fid, '# Domain: rectangle [%g, %g] x [%g, %g] (Chebfun spectral collocation)\n', a, b, c, d);
    fprintf(fid, '# Computed %s\n', datestr(now, 'yyyy-mm-dd HH:MM:SS'));
    fprintf(fid, '# N = %d Chebyshev points per direction, dofs = %d\n', N, numel(eigs_chebfun));
    fclose(fid);

    n        = (1:numel(eigs_chebfun))';
    lambda_n = eigs_chebfun(:);
    writetable(table(n, lambda_n), csv_file, 'WriteMode', 'append', 'WriteVariableNames', true);

    fprintf('Wrote %d eigenvalues to %s\n', numel(eigs_chebfun), csv_file);
end
