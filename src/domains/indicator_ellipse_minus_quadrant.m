function res = indicator_ellipse_minus_quadrant(x, y, a, b, c, e)

    tol = 1e-14;


    res1 = double(x > a & x < b & y > -(1/2)*sqrt(4-x.^2)+tol & y < e);

    res2 = double(x > a & x < c & y > e & y < (1/2)*sqrt(4-x.^2)-tol);

    % Adding interior boundary

    resb = double(x > a & x < b & y == e);

    res = res1 | res2 | resb;

end