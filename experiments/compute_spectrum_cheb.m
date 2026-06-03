function compute_spectrum_cheb(name, Nvals)
%COMPUTE_SPECTRUM_CHEB Run the Chebfun rectangle eigenvalue experiments.
%
%   compute_spectrum_cheb() computes the Dirichlet-Laplacian spectrum for every
%   domain in DOMAIN_CATALOG_CHEB at each resolution in a default list of
%   Chebyshev orders, writing one CSV per (domain, N) into results/eigenvalues/cheb/.
%
%   compute_spectrum_cheb(name) restricts to the single domain "name".
%   compute_spectrum_cheb(name, Nvals) overrides the list of Chebyshev orders N.
%   Pass "" or [] for name to keep all domains while choosing Nvals.
%
%   Examples:
%     compute_spectrum_cheb();                 % all domains, default Nvals
%     compute_spectrum_cheb("square");         % one domain, default Nvals
%     compute_spectrum_cheb("", [16 24 32]);   % all domains, chosen Nvals
%
%   Each CSV (columns n, lambda_n) carries a header recording the domain, box,
%   N, and dofs. Requires Chebfun; see startup.m.
%
%   See also DOMAIN_CATALOG_CHEB, CHEBFUN_LAPLACE_SPECTRUM.

    if nargin < 1
        name = '';
    end
    if nargin < 2 || isempty(Nvals)
        Nvals = [10 20 30 40];
    end

    % This file lives in experiments/; put the project sources and this folder
    % (for domain_catalog_cheb) on the path.
    experiments_dir = fileparts(mfilename('fullpath'));
    project_root    = fileparts(experiments_dir);
    run(fullfile(project_root, 'startup.m'));
    addpath(experiments_dir);

    out_dir = fullfile(project_root, 'results', 'eigenvalues', 'cheb');
    if ~exist(out_dir, 'dir')
        mkdir(out_dir);
    end

    cases = domain_catalog_cheb(name);
    for i = 1:numel(cases)
        c = cases(i);
        for N = Nvals
            fprintf('=== %s : N = %d ===\n', c.name, N);
            [evals, info] = chebfun_laplace_spectrum(c.box, N);
            csv_file = fullfile(out_dir, sprintf('%s_N%d-eigenvalues.csv', c.name, N));
            write_cheb_csv(csv_file, evals, c, info);
        end
    end
end

function write_cheb_csv(csv_file, evals, c, info)
    box = c.box;
    fid = fopen(csv_file, 'w');
    if fid == -1
        error('compute_spectrum_cheb:cannotOpen', 'Could not open %s for writing.', csv_file);
    end
    fprintf(fid, '# Domain: %s [%g, %g] x [%g, %g] (Chebfun spectral collocation)\n', ...
        c.name, box(1), box(2), box(3), box(4));
    fprintf(fid, '# Computed %s\n', datestr(now, 'yyyy-mm-dd HH:MM:SS'));
    fprintf(fid, '# N = %d Chebyshev points per direction, dofs = %d\n', info.N, info.dofs);
    if isfield(info, 'time')
        fprintf(fid, '# Computation time: %.3f s\n', info.time);
    end
    fclose(fid);

    n        = (1:numel(evals))';
    lambda_n = evals(:);
    writetable(table(n, lambda_n), csv_file, 'WriteMode', 'append', 'WriteVariableNames', true);
end
