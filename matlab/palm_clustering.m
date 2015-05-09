% Performs clustering based on PALM approach over columns of @A.
% Does exactly @iters iterations, searches for @k number of centers,
% assumes the data matrix @A is in R^(n x m).
% Returns @X, where for each iteration t, X(:,:,t) are the calculated
% centers, W(:,:,t) are the coefficients, and Psi(t) is the value of Psi.
function [X,W,Psi] = palm_clustering(A,n,m,k,iters)
    X = zeros(n,k,(iters+1));
    W = zeros(k,m,(iters+1));
    ones_vec = ones(m,1);
    Psi = zeros(1,iters);
    
    % X init
    %rand_vector_ind = randperm(m);
    %X(:,:,1) = A(:, rand_vector_ind(1:k));
    X(:,:,1) = clustering_init(A,n,m,k);
    
    % W init
    for i = 1:m
        W(:,i,1) = projection_onto_simplex(rand(k,1));
    end

    for t = 1:iters
        
        % W update
        alpha = min(W(:,:,t)*ones_vec);
        if alpha < 0.0000001
            disp(['alpha close to zero: ', num2str(alpha)]);
        end
        for i = 1:m
            d = distance_like(X(:,:,t), A(:,i), k);
            W(:,i,t+1) = projection_onto_simplex(W(:,i,t) - d/alpha);
        end
        
        % X update
        for l = 1:k
            W_l_row = W(l,:,t+1);
            X(:,l,t+1) = (A*W_l_row')/sum(W_l_row);
        end
        
        % Psi computations
        for i = 1:m
            Psi(t) = Psi(t) + distance_like(X(:,:,t+1), A(:,i), k)'*W(:,i,t+1);
        end
    end
    
end

