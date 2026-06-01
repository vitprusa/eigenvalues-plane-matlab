function compute_spectrum_fem_solvepdeeig(name)
%COMPUTE_SPECTRUM_FEM_SOLVEPDEEIG Run the FEM solvepdeeig eigenvalue experiments.
%
%   compute_spectrum_fem_solvepdeeig() computes the Dirichlet-Laplacian
%   spectrum for every domain in DOMAIN_CATALOG_FEM with the high-level
%   solvepdeeig solver over the range [0, 200]. It writes one CSV per domain
%   into results/fem/ as <domain>_fem_solvepdeeig-eigenvalues.csv (columns n,
%   lambda_n) with a metadata header.
%
%   compute_spectrum_fem_solvepdeeig(name) restricts to the single domain
%   "name".
%
%   Requires the PDE Toolbox.
%
%   See also DOMAIN_CATALOG_FEM, FEM_LAPLACE_SPECTRUM, COMPUTE_SPECTRUM_FEM_EIG.

    if nargin < 1
        name = '';
    end

    % This file lives in experiments/; put the project sources and this folder
    % (for domain_catalog_fem) on the path.
    experiments_dir = fileparts(mfilename('fullpath'));
    project_root    = fileparts(experiments_dir);
    run(fullfile(project_root, 'startup.m'));
    addpath(experiments_dir);

    out_dir = fullfile(project_root, 'results', 'fem');
    if ~exist(out_dir, 'dir')
        mkdir(out_dir);
    end

    cases = domain_catalog_fem(name);
    for i = 1:numel(cases)
        c = cases(i);
        fprintf('=== %s : FEM solvepdeeig ===\n', c.name);
        [evals, info] = fem_laplace_spectrum(c, "solvepdeeig");
        csv_file = fullfile(out_dir, sprintf('%s_fem_solvepdeeig-eigenvalues.csv', c.name));
        write_fem_csv(csv_file, evals, c, info);
    end
end
