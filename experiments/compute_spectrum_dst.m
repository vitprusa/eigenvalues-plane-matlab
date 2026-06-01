function compute_spectrum_dst(name, modes)
%COMPUTE_SPECTRUM_DST Run the DST-Laplacian eigenvalue experiments.
%
%   compute_spectrum_dst() computes both the full and the partial spectrum for
%   every domain in DOMAIN_CATALOG_DST and writes one CSV per (domain, mode) into
%   results/dst/.
%
%   compute_spectrum_dst(name) restricts the run to the single domain "name".
%   compute_spectrum_dst(name, modes) further restricts the modes, where modes
%   is "full", "partial", or ["full" "partial"]. Pass "" or [] for name to
%   keep all domains while selecting modes.
%
%   Examples:
%     compute_spectrum_dst();                        % everything
%     compute_spectrum_dst("L_shaped");              % one domain, both modes
%     compute_spectrum_dst("L_shaped", "partial");   % one domain, one mode
%     compute_spectrum_dst("", "full");              % all domains, full only
%
%   See also DOMAIN_CATALOG_DST, DST_LAPLACE_SPECTRUM.

    if nargin < 1
        name = '';
    end
    if nargin < 2 || isempty(modes)
        modes = ["full" "partial"];
    end
    modes = string(modes);

    % This file lives in experiments/; put the project sources and this
    % folder (for domain_catalog_dst) on the path.
    experiments_dir = fileparts(mfilename('fullpath'));
    project_root    = fileparts(experiments_dir);
    run(fullfile(project_root, 'startup.m'));
    addpath(experiments_dir);

    out_dir = fullfile(project_root, 'results', 'dst');

    % Partial-spectrum solver options, shared by all domains.
    opts = struct('k', 10, 'tolerance', 1e-10, ...
                  'subspace_dim', 100, 'max_iterations', 300);

    cases = domain_catalog_dst(name);
    for i = 1:numel(cases)
        for mode = modes
            fprintf('=== %s : %s spectrum ===\n', cases(i).name, mode);
            dst_laplace_spectrum(cases(i), mode, out_dir, opts);
        end
    end
end
