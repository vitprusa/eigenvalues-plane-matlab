function d2valsvec = dst_d2_chunk(valsvec, h)
%DST_D2_CHUNK Second derivative via discrete sine transform (DST).
%
%   d2valsvec = dst_d2_chunk(valsvec, h)
%
%   Computes the second derivative on a contiguous block of interior grid
%   points using the discrete sine transform (DST). Homogeneous zero
%   Dirichlet boundary conditions are assumed at the endpoints.
%
%   Inputs:
%     valsvec   - vector of function values (no NaNs) on a uniform grid
%                 with spacing h.
%     h         - grid spacing (scalar)
%
%   Output:
%     d2valsvec - vector of second derivatives at the same grid points.

    N = length(valsvec);
    L = (N + 1)*h;
    s = (pi/L)^2;
    kvec = -((1:N).^2)';
    d2valsvec = s * idst(kvec .* dst(valsvec(:)));

end
