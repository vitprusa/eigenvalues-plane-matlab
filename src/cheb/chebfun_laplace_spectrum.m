function [evals, info] = chebfun_laplace_spectrum(box, N)
%CHEBFUN_LAPLACE_SPECTRUM Dirichlet-Laplacian spectrum on a rectangle via Chebfun.
%
%   [evals, info] = chebfun_laplace_spectrum(box, N)
%
%   Computes the eigenvalues of the Dirichlet Laplacian on the rectangle
%   box = [a b c d] = [x_min x_max y_min y_max] by Chebyshev spectral
%   collocation: second-derivative matrices from Chebfun's diffmat in each
%   direction, interior points only (homogeneous Dirichlet), combined into the
%   2-D Laplacian with a Kronecker sum, then a dense eig.
%
%   Inputs:
%     box - [a b c d] rectangle extent.
%     N   - number of Chebyshev points per direction (default 20).
%
%   Outputs:
%     evals - (N-2)^2-by-1 eigenvalues of -Laplacian, ascending. Spectral
%             collocation resolves the low modes accurately; the high modes
%             are spurious.
%     info  - struct with fields N, box, dofs.
%
%   Requires Chebfun (diffmat) on the path; see startup.m.
%
%   See also COMPUTE_CHEB_SPECTRA, CHEB_DOMAIN_CATALOG.

    if nargin < 2 || isempty(N)
        N = 20;
    end
    a = box(1); b = box(2); c = box(3); d = box(4);

    % Second-derivative matrices on N Chebyshev points in each direction.
    D2x = diffmat(N, 2, [a b]);
    D2y = diffmat(N, 2, [c d]);

    % Drop the first/last points for homogeneous Dirichlet conditions.
    D2x_int = D2x(2:end-1, 2:end-1);
    D2y_int = D2y(2:end-1, 2:end-1);

    % 2-D Laplacian via Kronecker sum.
    I = speye(N - 2);
    L = kron(I, D2x_int) + kron(D2y_int, I);

    % Eigenvalues of -Laplacian (positive), ascending. Real part as a hedge
    % against spurious imaginary parts.
    evals = sort(real(eig(full(-L))), 'ascend');

    info = struct('N', N, 'box', box, 'dofs', (N - 2)^2);
end
