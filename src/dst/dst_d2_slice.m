function d2valsvec = dst_d2_slice(valsvec, h)
%DST_D2_SLICE Compute second derivative on a 1-D slice using the DST.
%
%   d2valsvec = dst_d2_slice(valsvec, h)
%
%   Computes the second derivative on a 1-D slice using the discrete sine
%   transform (DST). NaN values in the input mark points outside the
%   domain. Each contiguous block of non-NaN values is treated as an
%   independent segment with homogeneous Dirichlet boundary conditions at
%   its endpoints.
%
%   Inputs:
%     valsvec   - vector of function values on a uniform grid with spacing
%                 h. Points outside the domain are marked with NaN.
%     h         - grid spacing (scalar)
%
%   Output:
%     d2valsvec - vector of the same size as valsvec. Non-NaN entries
%                 contain the second derivative; NaN entries are preserved.

    mask = ~isnan(valsvec);
    d2valsvec = NaN(size(valsvec));

    if ~any(mask)
        return
    end

    % Find contiguous blocks of non-NaN values.
    %
    % The mask is padded with false on both sides and diff is applied.
    % Transitions false->true (+1) mark block starts, true->false (-1)
    % mark one past block ends.
    %
    % Example:
    %   valsvec = [2.1, 0.5, NaN, NaN, 3.7, 1.2, 0.8]
    %   mask    = [  1,   1,   0,   0,   1,   1,   1 ]
    %   padded  = [0, 1, 1, 0, 0, 1, 1, 1, 0]
    %   edges   = [1, 0, -1, 0, 1, 0, 0, -1]
    %   block_starts = [1, 5],  block_ends = [2, 7]
    %   => two blocks: indices 1:2 and 5:7
    edges = diff([false, mask(:)', false]);
    block_starts = find(edges == 1);
    block_ends = find(edges == -1) - 1;

    % Apply DST-based second derivative to individual blocks
    for k = 1:length(block_starts)
        idx = block_starts(k):block_ends(k);
        d2valsvec(idx) = dst_d2_chunk(valsvec(idx), h);
    end

end
