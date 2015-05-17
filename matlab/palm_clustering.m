% Performs clustering based on PALM approach over columns of @A.
% Does at most @max_iters iterations, searches for @k number of centers,
% assumes the data matrix @A is in R^(n x m).
% Returns @X, where for each iteration t, X(:,:,t) are the calculated
% centers, W(:,:,t) are the coefficients, Psi(t) is the value of Psi, I(t)
% are the clusters, and t is the number of iterations actual done.
function [X,I,t] = palm_clustering(A,n,m,k,max_iters,tol)
    X = zeros(n,k,(max_iters+1));
    W = zeros(k,m,(max_iters+1));
    I = zeros(m,(max_iters+1));
    Psi = zeros(1,max_iters);
    ones_vec = ones(m,1);
    u = 0.000001; % alpha multiplier
    
    % X init & init clustering
    X(:,:,1) = clustering_init(A,n,m,k);
    [~,CIDX] = clustering_distance(X(:,:,1), A, m, k);
    I(:,1) = CIDX';
    
    % W init
    for i = 1:m
        W(:,i,1) = projection_onto_simplex(rand(k,1));
    end

    for t = 1:max_iters
        
        % W update
        alpha = u*min(W(:,:,t)*ones_vec);
        if alpha < 0.0000001
            disp(['alpha close to zero: ', num2str(alpha)]);
        end
        for i = 1:m
            d = distance_like(X(:,:,t), A(:,i), k);
            W(:,i,t+1) = projection_onto_simplex(W(:,i,t) - 10*d/alpha);
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
        
        % clustering update
        [~,CIDX] = clustering_distance(X(:,:,t+1), A, m, k);
        I(:,t+1) = CIDX';
        
        if ((sum(I(:,t+1) == I(:,t)) == m) && t>1 && (Psi(t-1)-Psi(t))<tol)
            break;
        end
    end
    
end

