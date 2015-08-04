function [res] = H_l(w_l,x_l,A,m)
    res = 0;
    
    for i = 1:m
        res = res + w_l(i)*norm(x_l - A(:,i));
    end
end