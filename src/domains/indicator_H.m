function res = indicator_H(x, y)
%INDICATOR_H Indicator function for the letter H composed of seven unit squares.
%
%   res = indicator_H(x, y)
%
%   The letter H is composed of seven 1x1 squares:
%     - left column:  3 squares at x in (-1, 0), y in (-2, 1)
%     - center bar:   1 square  at x in ( 0, 1), y in (-1, 0)
%     - right column: 3 squares at x in ( 1, 2), y in (-2, 1)
%
%   The lower-left corner of the letter is at [-1, -2].
%
%   Inputs:
%     x, y - scalar or array coordinates (same size)
%
%   Output:
%     res  - logical array, true inside the domain

    res_left  = double(x > -1 & x < 0 & y > -2 & y < 1);
    res_bar   = double(x >  0 & x < 1 & y > -1 & y < 0);
    res_right = double(x >  1 & x < 2 & y > -2 & y < 1);

    % Interior boundaries between adjacent sub-regions
    res_b = double((x == 0 & y > -1 & y < 0) ...
                 | (x == 1 & y > -1 & y < 0));

    res = res_left | res_bar | res_right | res_b;

end
