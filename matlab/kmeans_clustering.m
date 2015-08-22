% Performs clustering based on KMEANS approach over columns of @A.
% Does at most @max_iters iterations, searches for @k number of centers,
% assumes the data matrix @A is in R^(n x m).
% Returns @X, where for each iteration t, X(:,:,t) are the calculated
% centers, and I(t) are the clusters indeces.
function [X,I,t,Phi] = kmeans_clustering(A,n,m,k,max_iters,tol,X_0)
    setenv('distance', 'sq-E-norm');
    X = zeros(n,k,(max_iters+1));
    I = zeros(m,(max_iters+1));
    Phi = zeros(1,(max_iters+1));
    ones_vec = ones(m,1);
    
    % X init & init clustering
    X(:,:,1) = X_0;
    [D,CIDX] = clustering_distance(X(:,:,1), A, m, k);
    I(:,1) = CIDX';
    Phi(1) = D*ones_vec;

    for t = 1:max_iters
        
        % X update
        for l = 1:k
            CIDX_labeled_l = (I(:,t)==l);
            if (sum(CIDX_labeled_l) == 0)
                disp(['no points in center: ', num2str(l)]);
            end
            X(:,l,t+1) = (A*CIDX_labeled_l)/sum(CIDX_labeled_l);
        end
        
       % I update
        [D,CIDX] = clustering_distance(X(:,:,t+1), A, m, k);
        I(:,t+1) = CIDX';
        Phi(t+1) = D*ones_vec;
        
        if ((sum(I(:,t+1) == I(:,t)) == m) && abs(Phi(t+1)-Phi(t))<tol)
            break;
        end
    end
    
    X = X(:,:,[1:t+1]);
    I = I(:,[1:t+1]);
    Phi = Phi(:,[1:t+1]);
end

