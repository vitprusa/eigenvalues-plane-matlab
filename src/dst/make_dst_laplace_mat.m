function [L, info] = make_dst_laplace_mat(x_range, y_range, M, indicator_function)
%MAKE_DST_LAPLACE_MAT Assemble the masked DST-based discrete Laplace matrix.
%
%   [L, info] = make_dst_laplace_mat(x_range, y_range, M, indicator_function)
%
%   Wrapper around MAKE_DST_LAPLACE_OP that materialises the operator as a
%   dense matrix. The i-th column of L is the operator applied to the i-th
%   column of the identity in the degrees-of-freedom space, so that
%   L * v == dst_laplace_op(v) for any DOF vector v.
%
%   The signature matches MAKE_DST_LAPLACE_OP; see that function for a
%   description of the inputs and of the returned info struct.
%
%   Outputs:
%     L    - info.dofs-by-info.dofs matrix of the discrete Laplacian on the
%            degrees of freedom.
%     info - grid/domain struct returned by MAKE_DST_LAPLACE_OP.
%
%   See also MAKE_DST_LAPLACE_OP, DST_LAPLACE_GRID.

    [dst_laplace_op, info] = make_dst_laplace_op(x_range, y_range, M, indicator_function);

    dofs = info.dofs;

    % Assemble matrix whose i-th column is the operator applied to the i-th
    % unit vector.
    L = zeros(dofs, dofs);
    parfor i = 1:dofs
        e = zeros(dofs, 1);
        e(i) = 1;
        L(:, i) = dst_laplace_op(e);
    end

end
