function cases = domain_catalog_dst(name)
%DOMAIN_CATALOG_DST Registry of domains for the DST-Laplacian experiments.
%
%   cases = DOMAIN_CATALOG_DST() returns a struct array describing every domain.
%   cases = DOMAIN_CATALOG_DST(name) returns only the entry whose name matches,
%   erroring if there is no such domain.
%
%   Each entry has the fields:
%     name      - char identifier, used in the output file names.
%     box       - bounding box [a b c d] = [x_min x_max y_min y_max].
%     phi       - indicator function handle phi(x, y), positive inside the
%                 domain. The domain must be embedded in the bounding box.
%     M_full    - resolution for the full-spectrum (dense eig) computation.
%     M_partial - resolution for the partial-spectrum (iterative eigs) one.
%
%   The M values are chosen so that the domain corners and edges fall on grid
%   lines; the relevant alignment condition is noted per row below.
%
%   To add a domain, add one make_case row here -- no new scripts are needed.
%
%   See also DST_LAPLACE_SPECTRUM, COMPUTE_SPECTRUM_DST.

    % Start with an empty struct array carrying the right fields.
    cases = repmat(make_case('', [0 0 0 0], @(x, y) x, 0, 0), 0, 1);

    % name                         box [a b c d]   indicator                                                    M_full  M_partial
    cases(end+1) = make_case('square',                 [0 pi 0 pi],   @(x, y) indicator_rectangle(x, y, 0, pi, 0, pi),               50, 200);  % fills the bounding box; analytic spectrum m^2 + n^2
    cases(end+1) = make_case('rectangle',              [0 2*pi 0 pi], @(x, y) indicator_rectangle(x, y, 0, 2*pi, 0, pi),            49, 199);  % fills the box; M+1 even so top edge y=pi on grid; analytic m^2/4 + n^2
    cases(end+1) = make_case('isosceles_triangle',     [0 pi 0 pi],   @(x, y) indicator_isosceles_triangle(x, y, 0, pi, 0),          50, 200);
    cases(end+1) = make_case('small_rectangle',        [0 pi 0 pi], @(x, y) indicator_rectangle(x, y, 0, pi/2, 0, pi/4),           47, 199);  % M+1 multiple of 4: edges pi/2, pi/4 on grid
    cases(end+1) = make_case('L_shaped',               [-1 1 -1 1], @(x, y) indicator_L_shaped(x, y, -1, 0, 1, -1, 0, 1),          49, 299);  % M+1 even: re-entrant corner at origin on grid
    cases(end+1) = make_case('ellipse_minus_quadrant', [-2 2 -1 1], @(x, y) indicator_ellipse_minus_quadrant(x, y, -2, 0, 2, 0),   47, 199);  % M+1 divisible by 4: x=0 and y=0 cut edges on grid
    cases(end+1) = make_case('H_shaped',               [-1 2 -2 1], @(x, y) indicator_H_shaped(x, y),                              50, 200);
    cases(end+1) = make_case('gww1',                   [-3 3 -3 3], @(x, y) indicator_gww1(x, y, -3, -1, 1, 3, -1, 1),             50, 200);  % isospectral with gww2
    cases(end+1) = make_case('gww2',                   [-3 3 -3 3], @(x, y) indicator_gww2(x, y, -3, -1, 1, 3, -1, 1, 3),          50, 200);  % isospectral with gww1

    if nargin >= 1 && ~isempty(name)
        match = strcmp({cases.name}, name);
        if ~any(match)
            error('domain_catalog_dst:unknownDomain', ...
                'Unknown domain "%s". Known domains: %s.', name, strjoin({cases.name}, ', '));
        end
        cases = cases(match);
    end
end

function c = make_case(name, box, phi, M_full, M_partial)
    c = struct('name', name, 'box', box, 'phi', phi, ...
               'M_full', M_full, 'M_partial', M_partial);
end
