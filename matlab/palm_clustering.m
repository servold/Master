% Performs clustering based on PALM approach over columns of @A.
% Does at most @max_iters iterations, searches for @k number of centers,
% assumes the data matrix @A is in R^(n x m).
% Returns @X, where for each iteration t, X(:,:,t) are the calculated
% centers, W(:,:,t) are the coefficients, Psi(t) is the value of Psi, I(t)
% are the clusters, and t is the number of iterations actual done.
function [X,W,I,I_W,t,Psi] = palm_clustering(A,n,m,k,max_iters,tol,x_0,w_0,alpha_update_f)
    setenv('distance', 'sq-E-norm');
    X = zeros(n,k,(max_iters+1));
    W = zeros(k,m,(max_iters+1));
    I = zeros(m,(max_iters+1));
    I_W = zeros(m,(max_iters+1));
    Psi = zeros(1,max_iters+1);
    ones_vec = ones(m,1);
    alpha0 = diam(A,m);
    
    % X init & init clustering
    X(:,:,1) = x_0;
    [~,CIDX] = clustering_distance(X(:,:,1), A, m, k);
    I(:,1) = CIDX';
    
    % W init
    W(:,:,1) = w_0;
    for j = 1:m
        [~,b] = max(w_0(:,j));
        I_W(j,1) = b;
    end
    
    % Psi computations
    for i = 1:m
        Psi(1) = Psi(1) + distance_like(X(:,:,1), A(:,i), k)'*W(:,i,1);
    end
    
    for t = 1:max_iters
        alpha = max(1e-8,alpha_update_f(alpha0,t));
        
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
            Psi(t+1) = Psi(t+1) + distance_like(X(:,:,t+1), A(:,i), k)'*W(:,i,t+1);
        end
        
        % clustering update
        [~,CIDX] = clustering_distance(X(:,:,t+1), A, m, k);
        I(:,t+1) = CIDX';
        
        for j = 1:m
            [~,b] = max(W(:,j,t+1));
            I_W(j,t+1) = b;
        end
        
        if ((Psi(t)-Psi(t+1))<tol)
            break;
        end
    end
    
    X = X(:,:,1:t+1);
    W = W(:,:,1:t+1);
    I = I(:,1:t+1);
    I_W = I_W(:,1:t+1);
    Psi = Psi(:,1:t+1);
end

