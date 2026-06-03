function write_eigs_csv(csv_file, lambda, c, mode, M, info, opts)
%WRITE_EIGS_CSV Write Dirichlet-Laplacian eigenvalues to CSV with a header.
%
%   write_eigs_csv(csv_file, lambda, c, mode, M, info, opts)
%
%   Writes a CSV file that begins with a block of "#"-prefixed comment lines
%   recording the parameters of the computation, followed by a table with
%   columns n and lambda_n.
%
%   Inputs:
%     csv_file - full path of the output file (its directory must exist).
%     lambda   - vector of eigenvalues of -Laplacian (positive), ordered from
%                the smallest.
%     c        - catalog entry with fields name, box, phi (see DOMAIN_CATALOG_DST).
%     mode     - "full" or "partial"; the eigs options are logged only for the
%                partial spectrum.
%     M        - grid resolution actually used.
%     info     - grid/domain struct; uses info.h and info.dofs.
%     opts     - partial-spectrum options with fields k, tolerance,
%                subspace_dim, max_iterations. Ignored when mode is "full".
%
%   See also DST_LAPLACE_SPECTRUM, DOMAIN_CATALOG_DST.

    mode = string(mode);
    box  = c.box;

    fid = fopen(csv_file, 'w');
    if fid == -1
        error('write_eigs_csv:cannotOpen', 'Could not open %s for writing.', csv_file);
    end

    fprintf(fid, '# Domain: %s (%s spectrum)\n', c.name, mode);
    fprintf(fid, '# Bounding box [a, b] x [c, d] = [%g, %g] x [%g, %g]\n', box(1), box(2), box(3), box(4));
    fprintf(fid, '# Resolution M = %d, grid spacing h = %g, dofs = %d\n', M, info.h, info.dofs);
    if isfield(info, 'time')
        fprintf(fid, '# Computation time: %.3f s\n', info.time);
    end
    fprintf(fid, '# Indicator function phi = %s\n', func2str(c.phi));
    if mode == "partial"
        fprintf(fid, '# eigs: k = %d, tolerance = %g, subspace_dim = %d, max_iterations = %d\n', ...
            opts.k, opts.tolerance, opts.subspace_dim, opts.max_iterations);
    end

    fclose(fid);

    % Append the eigenvalue table after the header comment block.
    n        = (1:numel(lambda))';
    lambda_n = lambda(:);
    writetable(table(n, lambda_n), csv_file, 'WriteMode', 'append', 'WriteVariableNames', true);
end
