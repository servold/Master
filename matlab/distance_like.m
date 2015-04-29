function [d] = distance_like(X,a)
    dim = size(X);
    k = dim(2);
    d = zeros(k,1);
    
    for l = 1:k
        d(l) = (norm(X(:,l) - a))^2;
    end
end

