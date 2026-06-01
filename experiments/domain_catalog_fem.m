function cases = domain_catalog_fem(name)
%DOMAIN_CATALOG_FEM Registry of domains for the FEM (PDE Toolbox) experiments.
%
%   cases = DOMAIN_CATALOG_FEM() returns a struct array of every domain.
%   cases = DOMAIN_CATALOG_FEM(name) returns only the matching entry, erroring
%   if there is no such domain.
%
%   Each entry has the fields:
%     name             - char identifier, used in the output file names.
%     gd, ns, sf       - decsg geometry description (geometry matrix, name
%                        matrix, set formula), passed as decsg(gd, sf, ns).
%     Hmax_eig         - target mesh size for the dense-eig workflow (coarser,
%                        because eig(K, M) is a dense solve).
%     Hmax_solvepdeeig - target mesh size for the solvepdeeig workflow (finer).
%
%   The Hmax values are tuning parameters; the eig workflow uses the coarser
%   one to keep the dense generalized eigenproblem affordable.
%
%   See also FEM_LAPLACE_SPECTRUM, COMPUTE_SPECTRUM_FEM_EIG,
%   COMPUTE_SPECTRUM_FEM_SOLVEPDEEIG.

    % Start with an empty struct array carrying the right fields.
    cases = repmat(make_case('', [], '', '', 0, 0), 0, 1);

    % name                  geometry (decsg gd, ns, sf)                                           Hmax_eig  Hmax_solvepdeeig
    cases(end+1) = make_case('square',             rect_gd(0, pi, 0, pi),    char('R1')',        'R1',     0.10, 0.05);
    cases(end+1) = make_case('rectangle',          rect_gd(0, 2*pi, 0, pi),  char('R1')',        'R1',     0.15, 0.07);
    % L-shape: [-1,1]^2 minus the upper-right unit square [0,1]x[0,1] (area 3).
    cases(end+1) = make_case('L_shaped',           [rect_gd(-1, 1, -1, 1), rect_gd(0, 1, 0, 1)], char('R1', 'R2')', 'R1-R2', 0.08, 0.05);
    % Isosceles right triangle (0,0)-(pi,0)-(pi,pi).
    cases(end+1) = make_case('isosceles_triangle', [2; 3; 0; pi; pi; 0; 0; pi], char('T1')',     'T1',     0.10, 0.05);

    if nargin >= 1 && ~isempty(name)
        match = strcmp({cases.name}, name);
        if ~any(match)
            error('domain_catalog_fem:unknownDomain', ...
                'Unknown domain "%s". Known domains: %s.', name, strjoin({cases.name}, ', '));
        end
        cases = cases(match);
    end
end

function c = make_case(name, gd, ns, sf, Hmax_eig, Hmax_solvepdeeig)
    c = struct('name', name, 'gd', gd, 'ns', ns, 'sf', sf, ...
               'Hmax_eig', Hmax_eig, 'Hmax_solvepdeeig', Hmax_solvepdeeig);
end

function gd = rect_gd(x1, x2, y1, y2)
%RECT_GD decsg geometry column for the rectangle [x1, x2] x [y1, y2].
    gd = [3; 4; x1; x2; x2; x1; y1; y1; y2; y2];
end
