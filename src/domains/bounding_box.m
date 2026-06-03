function [x_range, y_range] = bounding_box(a, b, c, d)
%BOUNDING_BOX Axis-aligned bounding box of a 2D domain as x/y ranges.
%
%   [x_range, y_range] = bounding_box(a, b, c, d)
%
%   Encapsulates the construction of the range vectors that the
%   MAKE_DST_LAPLACE_* factories expect, from the scalar corner coordinates
%   used in the scripts.
%
%   Inputs:
%     a, b - left and right x-coordinates of the bounding box (a < b).
%     c, d - bottom and top y-coordinates of the bounding box (c < d).
%
%   Outputs:
%     x_range - [a, b], the x-extent of the bounding box.
%     y_range - [c, d], the y-extent of the bounding box.
%
%   See also MAKE_DST_LAPLACE_OP, MAKE_DST_LAPLACE_MAT_BATCHED.

    x_range = [a, b];
    y_range = [c, d];

end
