function write_fd_csv(csv_file, evals, c, info)
%WRITE_FD_CSV Write finite-difference eigenvalues to CSV with a metadata header.
%
%   write_fd_csv(csv_file, evals, c, info)
%
%   Writes a CSV file that begins with a block of "#"-prefixed comment lines
%   recording the domain, bounding box, resolution, and grid spacing, followed
%   by a table with columns n and lambda_n.
%
%   Inputs:
%     csv_file - full path of the output file (its directory must exist).
%     evals    - vector of eigenvalues of -Laplacian, ascending.
%     c        - catalog entry (uses c.name and c.box).
%     info     - struct from FD_LAPLACE_SPECTRUM (M, h, dofs).
%
%   See also FD_LAPLACE_SPECTRUM, DOMAIN_CATALOG_FD.

    box = c.box;
    fid = fopen(csv_file, 'w');
    if fid == -1
        error('write_fd_csv:cannotOpen', 'Could not open %s for writing.', csv_file);
    end
    fprintf(fid, '# Domain: %s [%g, %g] x [%g, %g] (finite differences, 5-point stencil)\n', ...
        c.name, box(1), box(2), box(3), box(4));
    fprintf(fid, '# Computed %s\n', datestr(now, 'yyyy-mm-dd HH:MM:SS'));
    fprintf(fid, '# M = %d, grid spacing h = %g, dofs = %d\n', info.M, info.h, info.dofs);
    fclose(fid);

    n        = (1:numel(evals))';
    lambda_n = evals(:);
    writetable(table(n, lambda_n), csv_file, 'WriteMode', 'append', 'WriteVariableNames', true);
end
