
function res = indicator_Lshape(x, y, a, b, c, d, e, f)

    res1 = double(x > a & x < b & y > d & y < f);

    res2 = double(x >= b & x < c & y > d & y < e);

    % for i = 1:length(x)
    % 
    %     for j = 1:length(y)
    % 
    %         if (res1(i,j) > 0 | res2(i,j) > 0)
    % 
    %             res(i,j) = 1;
    % 
    %         else 
    % 
    %             res(i,j) = 0;
    % 
    %         end
    % 
    %     end
    % 
    % end

    res = res1 | res2;

end

    % if(x > a & x < b & y > d & y < f)
    % 
    %     s = 0;
    % 
    %     if(y < e)
    % 
    %         s = 1;
    % 
    %         while(s == 1)
    % 
    %             if(y > e)
    % 
    %                 s = 0;
    %                 res1 = double(x > a & x < b & y > d & y < f);
    % 
    %             end
    % 
    %         end
    %     end
    % 
    %     if(s == 1)
    % 
    %         res1 = double(x > a & x < b & y > d & y < e);
    % 
    %     else
    % 
    %         res1 = double(x > a & x < b & y > e & y < f);
    % 
    %     end
    % 
    % end
    % 
    % if(x >= b & x < c & y > d & y < f)
    % 
    %     s = 0;
    % 
    %     if(y < e)
    % 
    %         s = 1;
    % 
    %         while(s == 1)
    % 
    %             if(y > e)
    % 
    %                 s = 0;
    %                 res2 = double(x >= b & x < c & y > d & y < f);
    % 
    %             end
    % 
    %         end
    %     end
    % 
    %     if(s == 1)
    % 
    %         res2 = double(x >= b & x < c & y > d & y < e);
    % 
    %     else
    % 
    %         res2 = double(x >= b & x < c & y > e & y < f);
    % 
    %     end
    % end
    % 
    % res = [res1 res2];





