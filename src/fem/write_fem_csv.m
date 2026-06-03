function write_fem_csv(csv_file, evals, c, info)
%WRITE_FEM_CSV Write FEM eigenvalues to CSV with a metadata header.
%
%   write_fem_csv(csv_file, evals, c, info)
%
%   Writes a CSV file that begins with a block of "#"-prefixed comment lines
%   recording the domain, method, mesh size, and problem size, followed by a
%   table with columns n and lambda_n.
%
%   Inputs:
%     csv_file - full path of the output file (its directory must exist).
%     evals    - vector of eigenvalues of -Laplacian, ascending.
%     c        - catalog entry (uses c.name).
%     info     - struct from FEM_LAPLACE_SPECTRUM (method, Hmax, dofs, n_nodes).
%
%   See also FEM_LAPLACE_SPECTRUM, DOMAIN_CATALOG_FEM.

    fid = fopen(csv_file, 'w');
    if fid == -1
        error('write_fem_csv:cannotOpen', 'Could not open %s for writing.', csv_file);
    end
    fprintf(fid, '# Domain: %s (FEM, %s)\n', c.name, info.method);
    fprintf(fid, '# Computed %s\n', datestr(now, 'yyyy-mm-dd HH:MM:SS'));
    fprintf(fid, '# Hmax = %g, dofs = %d, mesh nodes = %d\n', info.Hmax, info.dofs, info.n_nodes);
    if isfield(info, 'time')
        fprintf(fid, '# Computation time: %.3f s\n', info.time);
    end
    fclose(fid);

    n        = (1:numel(evals))';
    lambda_n = evals(:);
    writetable(table(n, lambda_n), csv_file, 'WriteMode', 'append', 'WriteVariableNames', true);
end
