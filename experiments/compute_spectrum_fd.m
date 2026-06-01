function compute_spectrum_fd(name)
%COMPUTE_SPECTRUM_FD Run the finite-difference eigenvalue experiments.
%
%   compute_spectrum_fd() computes the Dirichlet-Laplacian spectrum for every
%   domain in DOMAIN_CATALOG_FD with the 5-point finite-difference stencil and
%   a dense eig, writing one CSV per domain into results/fd/ as
%   <domain>_fd-eigenvalues.csv (columns n, lambda_n) with a metadata header.
%
%   compute_spectrum_fd(name) restricts to the single domain "name".
%
%   See also DOMAIN_CATALOG_FD, FD_LAPLACE_SPECTRUM.

    if nargin < 1
        name = '';
    end

    % This file lives in experiments/; put the project sources and this folder
    % (for domain_catalog_fd) on the path.
    experiments_dir = fileparts(mfilename('fullpath'));
    project_root    = fileparts(experiments_dir);
    run(fullfile(project_root, 'startup.m'));
    addpath(experiments_dir);

    out_dir = fullfile(project_root, 'results', 'fd');
    if ~exist(out_dir, 'dir')
        mkdir(out_dir);
    end

    cases = domain_catalog_fd(name);
    for i = 1:numel(cases)
        c = cases(i);
        fprintf('=== %s : FD ===\n', c.name);
        [evals, info] = fd_laplace_spectrum(c);
        csv_file = fullfile(out_dir, sprintf('%s_fd-eigenvalues.csv', c.name));
        write_fd_csv(csv_file, evals, c, info);
    end
end
