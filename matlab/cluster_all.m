function [cluster] = cluster_all(X,A,m,k)
    cluster = zeros(1,m);
    
    for i = 1:m
       a = A(:,i);
       [~,cluster(i)] = min(distance_like(X,a,k));
    end
    
end

