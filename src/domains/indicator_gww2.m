function res = indicator_gww2(x, y, a, b, c, d, f, g, i)

    tol = 1e-15;

    res1 = ind_rect(x, y); 
    res2 = ind_iso_left_tri(x, y);
    res3 = ind_rev_iso_left_tri(x, y);

    % We have to include interior boundaries 
    res_b = double((x > a & x < b & y == g) | (x == c & y > f & y < g));
    

    res = res1 | res2 | res3 | res_b;
    
    
    function res_rect = ind_rect(x, y)

        res_rect = double(x > a & x < b & y > g & y < i);

    end

    function res_left_tri = ind_iso_left_tri(x, y)

        res_left_tri = double(x > a & x < c & y > (-x-2+tol) & y < g);

    end

    function res_reverse_left_tri = ind_rev_iso_left_tri(x, y)

        res_reverse_left_tri = double(x > c & x < d & y > f & y < (-x+2-tol));

    end

end







