function [mat] = concat_diff(a, b)
    length1 = size(a,2);
    length2 = size(b,2);
    l_diff = length1 - length2;
    
    if l_diff > 0
        b = [b, nan(1,l_diff)];
        mat = [a; b];
        
    elseif l_diff < 0
        a = [a, nan(size(a,1), abs(l_diff))];
        mat = [a; b];
       
    elseif l_diff == 0
        mat = [a; b];
    end
end