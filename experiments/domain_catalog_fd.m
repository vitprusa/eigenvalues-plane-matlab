function cases = domain_catalog_fd(name)
%DOMAIN_CATALOG_FD Registry of domains for the finite-difference experiments.
%
%   cases = DOMAIN_CATALOG_FD() returns a struct array of every domain.
%   cases = DOMAIN_CATALOG_FD(name) returns only the matching entry, erroring
%   if there is no such domain.
%
%   The domains match DOMAIN_CATALOG_FEM and DOMAIN_CATALOG_DST, but expressed
%   the way the finite-difference solver needs them: a rectangular bounding
%   box, an indicator function that masks the interior grid points, and a grid
%   resolution M -- the same representation the DST catalog uses.
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
    cases(end+1) = make_case('small_rectangle',        [0 pi 0 pi],   @(x, y) indicator_rectangle(x, y, 0, pi/2, 0, pi/4),         47);  % M+1 multiple of 4: edges pi/2, pi/4 on grid
    cases(end+1) = make_case('ellipse_minus_quadrant', [-2 2 -1 1],   @(x, y) indicator_ellipse_minus_quadrant(x, y, -2, 0, 2, 0), 47);  % M+1 divisible by 4: x=0 and y=0 cut edges on grid
    cases(end+1) = make_case('H',                      [-1 2 -2 1],   @(x, y) indicator_H(x, y),                                   50);
    cases(end+1) = make_case('gww1',                   [-3 3 -3 3],   @(x, y) indicator_gww1(x, y, -3, -1, 1, 3, -1, 1),           50);  % isospectral with gww2
    cases(end+1) = make_case('gww2',                   [-3 3 -3 3],   @(x, y) indicator_gww2(x, y, -3, -1, 1, 3, -1, 1, 3),        50);  % isospectral with gww1

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
