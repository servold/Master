function [X,W,d] = palm_clustering(A,n,m,k,iters)
    X = zeros(n,k);
    W = zeros(k,m);
    ones_v = ones(m,1);
    weights_v = ones(m,1);
    d = zeros(iters,1);
    
    rand_vector_ind = randperm(m);
    for l = 1:k
        X(:,l) = A(:, rand_vector_ind(l));
    end
    
    for i = 1:m
        W(:,i) = projection_onto_simplex(rand(k,1));
    end

    for t = 1:iters
        
        % W update
        alpha = min(W*ones_v);
        for i = 1:m
            d = distance_like(X, A(:,i));
            W(:,i) = projection_onto_simplex(W(:,i) - d/alpha);
        end
        
        % X update
        for l = 1:k
            X(:,l) = (A*W(l,:)')/sum(W(l,:));
        end
        
        d(t) = clustering_distance(X, A, weights_v, m, k);
    end
    
end

