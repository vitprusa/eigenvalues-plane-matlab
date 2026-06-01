function cases = domain_catalog_cheb(name)
%DOMAIN_CATALOG_CHEB Registry of rectangular domains for the Chebfun experiments.
%
%   cases = DOMAIN_CATALOG_CHEB() returns a struct array of every domain.
%   cases = DOMAIN_CATALOG_CHEB(name) returns only the matching entry, erroring
%   if there is no such domain.
%
%   Each entry has the fields:
%     name - char identifier, used in the output file names.
%     box  - bounding box [a b c d] = [x_min x_max y_min y_max].
%
%   The Chebfun method (diffmat + Kronecker sum) handles tensor-product
%   rectangles only. To add a domain, add one make_case row here.
%
%   See also COMPUTE_SPECTRUM_CHEB, CHEBFUN_LAPLACE_SPECTRUM.

    % Start with an empty struct array carrying the right fields.
    cases = repmat(make_case('', [0 0 0 0]), 0, 1);

    % name                 box [a b c d]
    cases(end+1) = make_case('square',    [0 pi   0 pi]);    % analytic n^2 + m^2
    cases(end+1) = make_case('rectangle', [0 2*pi 0 pi]);    % analytic m^2/4 + n^2

    if nargin >= 1 && ~isempty(name)
        match = strcmp({cases.name}, name);
        if ~any(match)
            error('domain_catalog_cheb:unknownDomain', ...
                'Unknown domain "%s". Known domains: %s.', name, strjoin({cases.name}, ', '));
        end
        cases = cases(match);
    end
end

function c = make_case(name, box)
    c = struct('name', name, 'box', box);
end
