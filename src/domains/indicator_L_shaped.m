
function res = indicator_L_shaped(x, y, a, b, c, d, e, f)

    res1 = double(x > a & x < b & y > d & y < f);

    res2 = double(x > b & x < c & y > d & y < e);

    resb = double(x == b & y > d & y < e);

    res = res1 | res2 | resb;

end

   
 





