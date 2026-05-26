function d2valsvec = dst_d2_slice(valsvec, h)
    % d2valsvec = dst_d2_slice(valsvec, h)
    %
    % Computes the second derivative on a 1-D slice using the discrete sine
    % transform (DST). NaN values in the input mark points outside the
    % domain; non-NaN values are the interior grid points on which the
    % derivative is computed (homogeneous Dirichlet boundary conditions are
    % assumed at the endpoints of each contiguous non-NaN segment).
    %
    % Inputs:
    %   valsvec  - vector of function values on a uniform grid with spacing
    %              h. Points outside the domain are marked with NaN.
    %   h        - grid spacing (scalar)
    %
    % Output:
    %   d2valsvec - vector of the same size as valsvec. Non-NaN entries
    %               contain the second derivative; NaN entries are preserved.

    mask = ~isnan(valsvec);
    d2valsvec = NaN(size(valsvec));

    if ~any(mask)
        return
    end

    interior = valsvec(mask);

    % Slice length, including the end ghost points
    L = (length(interior) + 1)*h;

    % Scaling for the second derivative
    s = (pi/L)^2;

    % Wavenumber vector for the second derivative
    kvec = -((1:length(interior)).^2)';

    % Compute the derivative at the interior points using DST
    d2valsvec(mask) = s * idst(kvec .* dst(interior(:)));
end
