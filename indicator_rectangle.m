% Returns 1 inside [a,b]x[c,d], 0 outside 
function res = indicator_rectangle(x, y, a, b, c, d)
    % x, y can be scalars or arrays (same size)
    res = double(x > a & x < b & y > c & y < d);
end
