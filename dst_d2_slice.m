function d2valsvec = dst_d2_slice(valsvec, h)
    % d2valsvec = dst_d2_slice(valsvec, h)
    %
    % Computes the second derivative of a 1-D slice using the discrete sine
    % transform (DST). The slice values are provided in valsvec and are
    % assumed to correspond to interior grid points (no ghost/end points).
    %
    % Inputs:
    %   valsvec  - column vector of function values on a uniform
    %              grid covering an interval of physical length L = (N+1)*h,
    %              where N = length(valsvec). The vector contains only the
    %              interior points (end/ghost points are omitted).
    %   h        - grid spacing (scalar)
    %
    % Output:
    %   d2valsvec - vector of second derivatives at the same interior grid
    %               points, computed with spectral accuracy using the DST
    %               (homogeneous Dirichlet boundary conditions are assumed).
    %
    % Notes:
    % - The routine rescales the spatial domain to the canonical interval
    %   [0, pi] before applying the DST-based second-derivative operator.
    % - The returned vector has the same orientation (row/column) as the
    %   input valsvec.
    %
    % Example:
    %   % Given interior values u_interior on a grid with spacing h:
    %   d2u = dst_d2_slice(u_interior, h);
    %

    % Slice length, including the end ghost points
    L = (length(valsvec) + 1)*h;
    
    % Scaling for the second derivative
    s = (pi/L)^2;
    
    % Wavenumber vector for the second derivative
    kvec = @(x) -((1:length(x)).^2);
    
    % Compute the derivative at the interior points using DST
    d2valsvec = s * idst(kvec(valsvec)' .* dst(valsvec));
end