function res = indicator_isosceles_triangle(x, y, a, b, c)

    res = double(x > a & x < b & y > c & y < x);

end


