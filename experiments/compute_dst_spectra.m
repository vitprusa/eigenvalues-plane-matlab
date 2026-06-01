function compute_dst_spectra(name, modes)
%COMPUTE_DST_SPECTRA Run the DST-Laplacian eigenvalue experiments.
%
%   compute_dst_spectra() computes both the full and the partial spectrum for
%   every domain in DOMAIN_CATALOG and writes one CSV per (domain, mode) into
%   results/dst/.
%
%   compute_dst_spectra(name) restricts the run to the single domain "name".
%   compute_dst_spectra(name, modes) further restricts the modes, where modes
%   is "full", "partial", or ["full" "partial"]. Pass "" or [] for name to
%   keep all domains while selecting modes.
%
%   Examples:
%     compute_dst_spectra();                        % everything
%     compute_dst_spectra("L_shaped");              % one domain, both modes
%     compute_dst_spectra("L_shaped", "partial");   % one domain, one mode
%     compute_dst_spectra("", "full");              % all domains, full only
%
%   See also DOMAIN_CATALOG, DST_LAPLACE_SPECTRUM.

    if nargin < 1
        name = '';
    end
    if nargin < 2 || isempty(modes)
        modes = ["full" "partial"];
    end
    modes = string(modes);

    % This file lives in experiments/; put the project sources and this
    % folder (for domain_catalog) on the path.
    experiments_dir = fileparts(mfilename('fullpath'));
    project_root    = fileparts(experiments_dir);
    run(fullfile(project_root, 'startup.m'));
    addpath(experiments_dir);

    out_dir = fullfile(project_root, 'results', 'dst');

    % Partial-spectrum solver options, shared by all domains.
    opts = struct('k', 10, 'tolerance', 1e-10, ...
                  'subspace_dim', 100, 'max_iterations', 300);

    cases = domain_catalog(name);
    for i = 1:numel(cases)
        for mode = modes
            fprintf('=== %s : %s spectrum ===\n', cases(i).name, mode);
            dst_laplace_spectrum(cases(i), mode, out_dir, opts);
        end
    end
end
