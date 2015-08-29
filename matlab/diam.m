function [d] = diam(A,m)
    d = 0;
    for i = 1:(m-1)
        for j = (i+1):m
            d = max(d, norm(A(:,i) - A(:,j)));
        end
    end
    
    sq_E_norm = strcmp(getenv('distance'),'sq-E-norm');
    if (sq_E_norm)
        d = d^2;
    end
end