% Computes the initial cluster centers. Starts with random center. Each
% iteration, adds the furthest data point from the centers picked so far.
function [X] = kmeans_pp_init(A,n,m,k)
    X = zeros(n,k);
    X(:,1) = A(:,randi(m,1,1));
    D = zeros(1,m);
    
    for l = 2:k
        partial_X = X(:,(1:(l-1)));
        for i = 1:m
            a = A(:,i);
            D(i) = min(distance_like(partial_X,a,(l-1)));
        end
        [~,j]= max(D);
        X(:,l) = A(:,j);
    end
end

