function [C] = cluster_all(X,A,m,k)
    C = zeros(m,1);
    
    for i = 1:m
       a = A(:,i);
       [~,cluster_index] = min(distance_like(X,a,k));
       C(i) = cluster_index;
    end
    
end

