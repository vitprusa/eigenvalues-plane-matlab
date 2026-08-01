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

    t_solve = tic;   % time the assembly and eigen-solve

    switch mode
        case "full"
            M = c.M_full;
            [L, info] = make_dst_laplace_mat_batched(x_range, y_range, M, c.phi);
            % Full spectrum via dense eig. The assembled matrix should be in
            % exact arithmetic symmetric, but in practice, for large assembled
            % matrices, the round-off errors in FFT lead to non-symmetric
            % matrix. The spurious non-symmetry that is a numerical artifact
            % then force EIG to use eigenvalues solver for non-symmetric
            % matrices. This is time consuming and gives eigenvalues with small
            % imaginary parts, which is again a numerical artifact. Previously
            % we have just applied real() to the so-obtained eigenvalues, but it
            % is better to symmetrise the assembled matrix beforehand, which
            % eliminates both the non-symmetric solver path and the spurious
            % imaginary parts in the computed eigenvalues. (Use
            % DST_LAPLACE_SYMMETRISE.) The impact of symmetrisation on
            % computation time and the computed eigenvalues is tested in
            % CHECK_LAPLACE_MATRIX_SYMMETRY, in test/laplace_matrix_symmetry.
            % lambda = sort(-real(eig(L)), 'ascend');
            lambda = sort(-eig(dst_laplace_symmetrise(L)), 'ascend');

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

    info.time = toc(t_solve);

    if ~exist(out_dir, 'dir')
        mkdir(out_dir);
    end
    csv_file = fullfile(out_dir, sprintf('%s_%s-eigenvalues.csv', c.name, mode));
    write_eigs_csv(csv_file, lambda, c, mode, M, info, opts);
end
