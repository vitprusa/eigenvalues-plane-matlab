function compute_spectrum_fem_eig(name)
%COMPUTE_SPECTRUM_FEM_EIG Run the FEM assembled-matrix eigenvalue experiments.
%
%   compute_spectrum_fem_eig() computes the Dirichlet-Laplacian spectrum for
%   every domain in DOMAIN_CATALOG_FEM by assembling the stiffness and mass
%   matrices and solving the dense generalized eigenproblem eig(K, M). It
%   writes one CSV per domain into results/fem/ as
%   <domain>_fem_eig-eigenvalues.csv (columns n, lambda_n) with a metadata
%   header.
%
%   compute_spectrum_fem_eig(name) restricts to the single domain "name".
%
%   Requires the PDE Toolbox.
%
%   See also DOMAIN_CATALOG_FEM, FEM_LAPLACE_SPECTRUM,
%   COMPUTE_SPECTRUM_FEM_SOLVEPDEEIG.

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
        fprintf('=== %s : FEM eig ===\n', c.name);
        [evals, info] = fem_laplace_spectrum(c, "eig");
        csv_file = fullfile(out_dir, sprintf('%s_fem_eig-eigenvalues.csv', c.name));
        write_fem_csv(csv_file, evals, c, info);
    end
end
