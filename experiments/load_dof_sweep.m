function results = load_dof_sweep(out_dir, domain, methods)
%LOAD_DOF_SWEEP Read the DOF-sweep CSVs for a domain into a results struct.
%
%   results = load_dof_sweep(out_dir, domain, methods) returns a struct with one
%   field per method; results.(method) is a cell array of
%   struct('evals', ..., 'dofs', ...), ordered by the numeric suffix k in the
%   file name <domain>_<method>_<k>-eigenvalues.csv.
%
%   See also READ_EIGS_CSV, PLOT_L_SHAPED_DOF_SWEEP, PLOT_RECTANGLE_DOF_SWEEP.

    results = struct();
    for mi = 1:numel(methods)
        m = methods{mi};
        files = dir(fullfile(out_dir, sprintf('%s_%s_*-eigenvalues.csv', domain, m)));
        ks = zeros(1, numel(files));
        for j = 1:numel(files)
            tok = regexp(files(j).name, '_(\d+)-eigenvalues\.csv$', 'tokens', 'once');
            ks(j) = str2double(tok{1});
        end
        [~, order] = sort(ks);
        files = files(order);

        cols = cell(1, numel(files));
        for j = 1:numel(files)
            [evals, dofs] = read_eigs_csv(fullfile(out_dir, files(j).name));
            cols{j} = struct('evals', evals, 'dofs', dofs);
        end
        results.(m) = cols;
    end
end
