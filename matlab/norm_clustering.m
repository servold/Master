function [X,I,t,Phi] = norm_clustering(A,n,m,k,max_iters,tol,x_0)
    setenv('distance', 'E-norm');
    X = zeros(n,k,(max_iters+1));
    W = zeros(k,m,(max_iters+1));
    I = zeros(m,(max_iters+1));
    Psi = zeros(1,max_iters);
    Phi = zeros(k,m,max_iters);
    ones_vec = ones(m,1);
    alpha0 = 100;
    
    [w0,x0] = SP_method(A,n,k,m);
    
    % x init & init clustering
    X(:,:,1) = x0;
    [~,CIDX] = clustering_distance(X(:,:,1), A, m, k);
    I(:,1) = CIDX';
    
    % W init
    W(:,:,1) = w0;
%     for i = 1:m
%         W(:,i,1) = rand(k,1);
%         W(:,i,1) = W(:,i,1)/sum(W(:,i,1));
%     end
    
    M = zeros(k,m);
    for l = 1:k
        H_l_w0_x0 = H_l(W(l,:,1), X(:,l,1), A,m);
        sum_w0_l = sum(W(l,:,1));
        for i = 1:m
            M(l,i) = (H_l(W(l,:,1), A(:,i), A,m) - H_l_w0_x0)/sum_w0_l;
        end
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
        
        norm_k_m = zeros(k,m);
        for l = 1:k
            for i = 1:m
                norm_k_m(l,i) = norm(X(:,l,t+1) - A(:,i));
            end
        end
        Phi(:,:,t) = norm_k_m - M;
        
        if ((sum(I(:,t+1) == I(:,t)) == m) && t>1 && (Psi(t-1)-Psi(t))<tol)
            break;
        end
    end
    
    X = X(:,:,[1:t+1]);
    W = W(:,:,[1:t+1]);
    I = I(:,[1:t+1]);
    Psi = Psi(:,[1:t]);
    Phi = Phi(:,:,[1:t]);
end
