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
    % Small rectangle [0, pi/2] x [0, pi/4].
    cases(end+1) = make_case('small_rectangle', rect_gd(0, pi/2, 0, pi/4), char('R1')', 'R1', 0.05, 0.03);
    % Ellipse (semi-axes 2, 1) minus the lower-right quadrant.
    cases(end+1) = make_case('ellipse_minus_quadrant', [ellipse_gd(0, 0, 2, 1), rect_gd(0, 2, -1, 0)], char('E1', 'R1')', 'E1-R1', 0.10, 0.06);
    % H: left, centre, right unit-square columns.
    cases(end+1) = make_case('H', [rect_gd(-1, 0, -2, 1), rect_gd(0, 1, -1, 0), rect_gd(1, 2, -2, 1)], char('R1', 'R2', 'R3')', 'R1+R2+R3', 0.08, 0.05);
    % GWW1 isospectral drum (two triangles + a rectangle).
    cases(end+1) = make_case('gww1', [tri_gd(-3, 1, -1, -1, -1, 3), rect_gd(-1, 3, -1, 1), tri_gd(1, -1, 1, -3, 3, -1)], char('TL', 'R1', 'TR')', 'TL+R1+TR', 0.15, 0.10);
    % GWW2 isospectral drum (a rectangle + two triangles).
    cases(end+1) = make_case('gww2', [rect_gd(-3, -1, 1, 3), tri_gd(-3, 1, 1, -3, 1, 1), tri_gd(1, 1, 1, -1, 3, -1)], char('R1', 'T2', 'T3')', 'R1+T2+T3', 0.15, 0.10);

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

function gd = tri_gd(x1, y1, x2, y2, x3, y3)
%TRI_GD decsg geometry column for the triangle with the given (CCW) vertices.
    gd = [2; 3; x1; x2; x3; y1; y2; y3; 0; 0];   % padded to length 10
end

function gd = ellipse_gd(xc, yc, ax, by)
%ELLIPSE_GD decsg geometry column for the axis-aligned ellipse at (xc, yc).
    gd = [4; xc; yc; ax; by; 0; 0; 0; 0; 0];     % padded to length 10
end
