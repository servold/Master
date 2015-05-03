% Computes a k-dimentional vector, each coordinate l is the distance 
% of the l column of @X from @a.
% Currently the applied distance is the squared euclidean norm.
function [d] = distance_like(X,a,k)
    d = zeros(k,1);
    
    for l = 1:k
        d(l) = (norm(X(:,l) - a))^2;
    end
end

