function [evals, info] = fd_laplace_spectrum(entry)
%FD_LAPLACE_SPECTRUM Dirichlet-Laplacian spectrum on a domain by finite differences.
%
%   [evals, info] = fd_laplace_spectrum(entry)
%
%   Computes the eigenvalues of the Dirichlet Laplacian on the domain
%   described by the catalog entry (see DOMAIN_CATALOG_FD) with the standard
%   5-point finite-difference stencil. The bounding box is sampled on a
%   uniform grid of spacing h = (b - a)/(M + 1) (the same h in both
%   directions); interior grid points are selected by the indicator function,
%   and the Laplacian couples each kept point to its in-domain neighbours, with
%   homogeneous Dirichlet conditions wherever a neighbour falls outside.
%
%   Inputs:
%     entry - struct with fields name, box = [a b c d], phi, M
%             (see DOMAIN_CATALOG_FD).
%
%   Outputs:
%     evals - dofs-by-1 eigenvalues of -Laplacian, ascending.
%     info  - struct with fields M, h, box, dofs.
%
%   See also DOMAIN_CATALOG_FD, COMPUTE_SPECTRUM_FD.

    t0 = tic;

    a = entry.box(1); b = entry.box(2); c = entry.box(3); d = entry.box(4);
    M = entry.M;
    h = (b - a) / (M + 1);

    x_vec = a : h : b;
    y_vec = c : h : d;
    [X, Y] = ndgrid(x_vec, y_vec);

    % Interior grid points selected by the indicator (the box boundary, where
    % the indicator is non-positive, is excluded).
    mask = entry.phi(X, Y) > 0;
    dofs = sum(mask, 'all');
    [nx, ny] = size(mask);

    % DOF index for each kept grid point (column-major, matching find).
    dof_map = zeros(nx, ny);
    dof_map(mask) = 1:dofs;
    [ii, jj] = find(mask);

    % Assemble the 5-point Laplacian as sparse triplets.
    rows = zeros(5*dofs, 1); cols = rows; vals = rows;
    cnt = 0;
    inv_h2 = 1 / h^2;
    for p = 1:dofs
        i = ii(p); j = jj(p);
        cnt = cnt + 1; rows(cnt) = p; cols(cnt) = p; vals(cnt) = -4 * inv_h2;
        if i > 1  && mask(i-1, j), cnt = cnt+1; rows(cnt)=p; cols(cnt)=dof_map(i-1, j); vals(cnt)=inv_h2; end
        if i < nx && mask(i+1, j), cnt = cnt+1; rows(cnt)=p; cols(cnt)=dof_map(i+1, j); vals(cnt)=inv_h2; end
        if j > 1  && mask(i, j-1), cnt = cnt+1; rows(cnt)=p; cols(cnt)=dof_map(i, j-1); vals(cnt)=inv_h2; end
        if j < ny && mask(i, j+1), cnt = cnt+1; rows(cnt)=p; cols(cnt)=dof_map(i, j+1); vals(cnt)=inv_h2; end
    end
    L = sparse(rows(1:cnt), cols(1:cnt), vals(1:cnt), dofs, dofs);

    % Eigenvalues of -Laplacian (positive), ascending. Real part as a hedge.
    evals = sort(real(eig(full(-L))), 'ascend');

    info = struct('M', M, 'h', h, 'box', entry.box, 'dofs', dofs, 'time', toc(t0));
end
