% Performs clustering based on PALM approach over columns of @A.
% Does at most @max_iters iterations, searches for @k number of centers,
% assumes the data matrix @A is in R^(n x m).
% Returns @X, where for each iteration t, X(:,:,t) are the calculated
% centers, W(:,:,t) are the coefficients, Psi(t) is the value of Psi, I(t)
% are the clusters, and t is the number of iterations actual done.
function [X,I,t] = palm_clustering(A,n,m,k,max_iters,tol,X_0,prefix)
    X = zeros(n,k,(max_iters+1));
    W = zeros(k,m,(max_iters+1));
    I = zeros(m,(max_iters+1));
    Psi = zeros(1,max_iters);
    ones_vec = ones(m,1);
    alpha = 100;
    
    % X init & init clustering
    X(:,:,1) = X_0;
    [~,CIDX] = clustering_distance(X(:,:,1), A, m, k);
    I(:,1) = CIDX';
    
    % plot x0
    fig = plot_clusters(A,CIDX,k,X_0);
%    saveas(fig,[prefix,'-alpha-',num2str(alpha),'-iter-',num2str(0),'.jpg']);
    
    % W init
    for i = 1:m
        W(:,i,1) = projection_onto_simplex(rand(k,1));
    end

    for t = 1:max_iters
        alpha = 100/t;
        
        % W update
        [v,j] = min(W(:,:,t)*ones_vec);
        if v < 0.0000001
            disp(['cluster ', num2str(j) ,' close to zero: ', num2str(v)]);
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
        
        % clustering update
        [~,CIDX] = clustering_distance(X(:,:,t+1), A, m, k);
        I(:,t+1) = CIDX';
        
        if (t <=10 || mod(t,10) == 1)
            fig = plot_clusters(A,CIDX,k,X(:,:,t+1));
%            saveas(fig,[prefix,'-alpha-',num2str(alpha),'-iter-',num2str(t),'.jpg']);
        end
        
        if ((sum(I(:,t+1) == I(:,t)) == m) && t>1 && (Psi(t-1)-Psi(t))<tol)
            break;
        end
    end
    
end

