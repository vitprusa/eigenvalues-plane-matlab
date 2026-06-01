function [lambda, info] = dst_laplace_spectrum(c, mode, out_dir, opts)
%DST_LAPLACE_SPECTRUM Compute and export a DST-Laplacian spectrum for one domain.
%
%   [lambda, info] = dst_laplace_spectrum(c, mode, out_dir, opts)
%
%   Computes the Dirichlet-Laplacian eigenvalues on the domain described by
%   the catalog entry c (see DOMAIN_CATALOG_DST) and writes them to
%   <out_dir>/<c.name>_<mode>-eigenvalues.csv with a metadata header.
%
%   The domain is described by an indicator function embedded in a rectangular
%   bounding box and discretised on a uniform grid. WARNING: the domain must
%   be embedded in the bounding box; this is not checked.
%
%   Inputs:
%     c       - struct with fields name, box = [a b c d], phi, M_full,
%               M_partial (see DOMAIN_CATALOG_DST).
%     mode    - "full"    : whole spectrum from a dense matrix
%                           (make_dst_laplace_mat_batched + eig), at M_full.
%               "partial" : leading opts.k eigenvalues from the matrix-free
%                           operator (make_dst_laplace_op_batched + eigs), at
%                           M_partial.
%     out_dir - directory for the output CSV; created if it does not exist.
%     opts    - partial-spectrum options with fields k, tolerance,
%               subspace_dim, max_iterations. Ignored when mode is "full".
%
%   Outputs:
%     lambda  - eigenvalues of -Laplacian, positive and ordered from the
%               smallest. The discrete operator approximates the (negative-
%               definite) Laplacian, so its eigenvalues are negated here.
%     info    - grid/domain struct returned by the operator factory.
%
%   See also DOMAIN_CATALOG_DST, WRITE_EIGS_CSV, MAKE_DST_LAPLACE_MAT_BATCHED,
%   MAKE_DST_LAPLACE_OP_BATCHED.

    if nargin < 4
        opts = struct();   % only needed for the partial spectrum
    end

    mode = string(mode);
    [x_range, y_range] = bounding_box(c.box(1), c.box(2), c.box(3), c.box(4));

    switch mode
        case "full"
            M = c.M_full;
            [L, info] = make_dst_laplace_mat_batched(x_range, y_range, M, c.phi);
            % Full spectrum via dense eig. Take the real part as a hedge
            % against spurious imaginary parts, then report -Laplacian.
            lambda = sort(-real(eig(L)), 'ascend');

        case "partial"
            M = c.M_partial;
            [Lop, info] = make_dst_laplace_op_batched(x_range, y_range, M, c.phi);
            % Leading eigenvalues via the iterative eigs. 'largestreal' picks
            % the eigenvalues of the Laplacian closest to zero (least negative),
            % i.e. the leading ones; the operator is symmetric.
            D = eigs(Lop, info.dofs, opts.k, 'largestreal', ...
                'Tolerance',           opts.tolerance, ...
                'MaxIterations',       opts.max_iterations, ...
                'SubspaceDimension',   opts.subspace_dim, ...
                'IsFunctionSymmetric', true, ...
                'Display',             1);
            lambda = sort(-real(D), 'ascend');

        otherwise
            error('dst_laplace_spectrum:badMode', ...
                'mode must be "full" or "partial", got "%s".', mode);
    end

    if ~exist(out_dir, 'dir')
        mkdir(out_dir);
    end
    csv_file = fullfile(out_dir, sprintf('%s_%s-eigenvalues.csv', c.name, mode));
    write_eigs_csv(csv_file, lambda, c, mode, M, info, opts);
end
