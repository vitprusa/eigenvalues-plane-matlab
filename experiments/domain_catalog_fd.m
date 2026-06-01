function cases = domain_catalog_fd(name)
%DOMAIN_CATALOG_FD Registry of domains for the finite-difference experiments.
%
%   cases = DOMAIN_CATALOG_FD() returns a struct array of every domain.
%   cases = DOMAIN_CATALOG_FD(name) returns only the matching entry, erroring
%   if there is no such domain.
%
%   The domains match DOMAIN_CATALOG_FEM (square, rectangle, L-shaped, and the
%   isosceles right triangle), but expressed the way the finite-difference
%   solver needs them: a rectangular bounding box, an indicator function that
%   masks the interior grid points, and a grid resolution M -- the same
%   representation the DST catalog uses.
%
%   Each entry has the fields:
%     name - char identifier, used in the output file names.
%     box  - bounding box [a b c d] = [x_min x_max y_min y_max].
%     phi  - indicator function handle phi(x, y), positive inside the domain.
%     M    - number of interior grid points along a full bounding-box slice;
%            the spacing is h = (b - a)/(M + 1), the same in both directions.
%
%   See also FD_LAPLACE_SPECTRUM, COMPUTE_SPECTRUM_FD, DOMAIN_CATALOG_FEM.

    % Start with an empty struct array carrying the right fields.
    cases = repmat(make_case('', [0 0 0 0], @(x, y) x, 0), 0, 1);

    % name                     box [a b c d]   indicator                                          M
    cases(end+1) = make_case('square',             [0 pi 0 pi],   @(x, y) indicator_rectangle(x, y, 0, pi, 0, pi),         50);
    cases(end+1) = make_case('rectangle',          [0 2*pi 0 pi], @(x, y) indicator_rectangle(x, y, 0, 2*pi, 0, pi),      49);  % M+1 even: top edge y=pi on grid
    cases(end+1) = make_case('L_shaped',           [-1 1 -1 1],   @(x, y) indicator_L_shaped(x, y, -1, 0, 1, -1, 0, 1),   49);  % M+1 even: re-entrant corner at origin on grid
    cases(end+1) = make_case('isosceles_triangle', [0 pi 0 pi],   @(x, y) indicator_isosceles_triangle(x, y, 0, pi, 0),    50);

    if nargin >= 1 && ~isempty(name)
        match = strcmp({cases.name}, name);
        if ~any(match)
            error('domain_catalog_fd:unknownDomain', ...
                'Unknown domain "%s". Known domains: %s.', name, strjoin({cases.name}, ', '));
        end
        cases = cases(match);
    end
end

function c = make_case(name, box, phi, M)
    c = struct('name', name, 'box', box, 'phi', phi, 'M', M);
end
