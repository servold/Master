% Computes a k-dimentional vector, each coordinate l is the distance 
% of the l column of @X from @a.
% Currently the applied distance is the squared euclidean norm.
function [d] = distance_like(X,a,k)
    sq_E_norm = strcmp(getenv('distance'),'sq-E-norm');
    d = zeros(k,1);
    
    for l = 1:k
        if (sq_E_norm)
            d(l) = (norm(X(:,l) - a))^2;
        else
            d(l) = (norm(X(:,l) - a));
        end
    end
end

