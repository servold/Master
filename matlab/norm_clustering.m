function [X,I,t] = norm_clustering(A,n,m,k,max_iters,tol,x_0)
    setenv('distance', 'E-norm');
    X = zeros(n,k,(max_iters+1));
    W = zeros(k,m,(max_iters+1));
    I = zeros(m,(max_iters+1));
    Psi = zeros(1,max_iters);
    ones_vec = ones(m,1);
    alpha0 = 100;
    
    % x init & init clustering
    X(:,:,1) = x_0;
    [~,CIDX] = clustering_distance(X(:,:,1), A, m, k);
    I(:,1) = CIDX';
    
    % W init
    for i = 1:m
        W(:,i,1) = rand(k,1);
        W(:,i,1) = W(:,i,1)/sum(W(:,i,1));
    end

    for t = 1:max_iters
        alpha = alpha0/t;
        
        % W update
        [v,j] = min(W(:,:,t)*ones_vec);
        if v < 0.0000001
            disp(['cluster ', num2str(j) ,' close to zero: ', num2str(v)]);
        end
        for i = 1:m
            d = distance_like(X(:,:,t), A(:,i), k);
            W(:,i,t+1) = projection_onto_simplex(W(:,i,t) - d/alpha);
        end
        
        % x update
        for l = 1:k
            x_l = X(:,l,t);
            w_l = W(l,:,t+1);
            u = zeros(m,1);
            
            for i = 1:m
                u(i) = w_l(i)/norm(x_l - A(:,i));
            end
            
            X(:,l,t+1) = A*u/sum(u);
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
    
    X = X(:,:,[1:t+1]);
    W = W(:,:,[1:t+1]);
    I = I(:,[1:t+1]);
    Psi = Psi(:,[1:t]);
end
