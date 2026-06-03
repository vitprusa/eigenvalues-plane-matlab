function res = indicator_gww1(x, y, a, b, c, d, f, g)

    tol = 1e-15;

    res1 = ind_iso_right_tri(x, y); 
    res2 = ind_iso_left_tri(x, y);
    res3 = ind_rect(x, y);
    res4 = ind_rev_iso_right_tri(x, y);

    % We have to include interior boundaries 
    res_b = double((x > a & x < b & y == g) | (x == b & y > f & y < g) ...
                    | (x > c & x < d & y == f));
    

    res = res1 | res2 | res3 | res4 | res_b;
    
    
    function res_rect = ind_rect(x, y)

        res_rect = double(x > b & x < d & y > f & y < g);

    end

    function res_right_tri = ind_iso_right_tri(x, y)

        res_right_tri = double(x > a & x < b & y > g & y < (x+4-tol));

    end

    function res_left_tri = ind_iso_left_tri(x, y)

        res_left_tri = double(x > a & x < b & y > (-x-2+tol) & y < g);

    end

    function res_reverse_right_tri = ind_rev_iso_right_tri(x, y)

        res_reverse_right_tri = double(x > c & x < d & y > (x-4+tol) & y < f);

    end

end







