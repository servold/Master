function [X,W,V,U,I,I_W,t,B_eps,H_eps,ch1,ch2] = eps_norm_clustering3(A,n,m,k,max_iters,tol,x_0,w_0,eps,alpha_update_f)
    setenv('distance', 'E-norm');
    X = zeros(n,k,(max_iters+1));
    W = zeros(k,m,(max_iters+1));
    V = zeros(k,m,(max_iters+1));
    U = zeros(k,m,(max_iters+1));
    I = zeros(m,(max_iters+1));
    I_W = zeros(m,(max_iters+1));
    B_eps = zeros(1,max_iters+1);
    H_eps = zeros(1,max_iters+1);
    ch1 = zeros(1,max_iters+1);
    ch2 = zeros(1,max_iters+1);
    ones_vec = ones(m,1);
    diameter = diam(A,m);
    kapa = sqrt(2)*max(diameter,eps);
    alpha0 = diameter;
    
    % x init & init clustering
    X(:,:,1) = x_0;
    [~,CIDX] = clustering_distance(X(:,:,1), A, m, k);
    I(:,1) = CIDX';
    
    % W init
    W(:,:,1) = w_0;
    [~,b] = max(w_0);
    I_W(:,1) = b';
    
    % V init
    V(:,:,1) = (1/sqrt(diameter^2 + eps^2))*ones(k,m);
    
    % B&H computations
    B_eps(1)=B_epsilon(A,m,k,eps,W(:,:,1),V(:,:,1),X(:,:,1));
    H_eps(1)=H_epsilon(A,m,k,eps,W(:,:,1),X(:,:,1));

    for t = 1:max_iters
        alpha = max(1e-8,alpha_update_f(alpha0,t));
        
        % W&V update
        [v,j] = min(W(:,:,t)*ones_vec);
        if v < 0.0000001 % sanity check
            disp(['cluster ', num2str(j) ,' close to zero: ', num2str(v)]);
        end
        for i = 1:m
            d_i = distance_like(X(:,:,t), A(:,i), k);
            d_i_eps = zeros(k,1);
            for l = 1:k
                V(l,i,t+1) = 1/sqrt(d_i(l)^2 + eps^2);
                d_i_eps(l) = sqrt(d_i(l)^2 + eps^2);
            end
            W(:,i,t+1) = projection_onto_simplex(W(:,i,t) - d_i_eps/alpha);
        end
        
        % X update
        for l = 1:k
            WV_l_row = W(l,:,t+1).*V(l,:,t+1);
            X(:,l,t+1) = (A*WV_l_row')/sum(WV_l_row);
        end
        
        % B&H computations
        B_eps(t+1)=B_epsilon(A,m,k,eps,W(:,:,t+1),V(:,:,t+1),X(:,:,t+1));
        H_eps(t+1)=H_epsilon(A,m,k,eps,W(:,:,t+1),X(:,:,t+1));
        ch1(t+1) = H_eps(t) - H_epsilon(A,m,k,eps,W(:,:,t+1),X(:,:,t));
        ch2(t+1) = H_epsilon(A,m,k,eps,W(:,:,t+1),X(:,:,t)) - H_eps(t+1);
        
        % clustering update
        [~,CIDX] = clustering_distance(X(:,:,t+1), A, m, k);
        I(:,t+1) = CIDX';
        [~,b] = max(W(:,:,t+1));
        I_W(:,t+1) = b';
        
        if ((H_eps(t)-H_eps(t+1))<tol)
            break;
        end
    end
    
    X = X(:,:,1:t+1);
    W = W(:,:,1:t+1);
    V = V(:,:,1:t+1);
    U = U(:,:,1:t+1);
    I = I(:,1:t+1);
    I_W = I_W(:,1:t+1);
    B_eps = B_eps(:,1:t+1);
    H_eps = H_eps(:,1:t+1);
end

function [res] = B_epsilon(A,m,k,eps,w,v,x)
    setenv('distance', 'E-norm');
    res = 0;
    for i = 1:m
        d_i = distance_like(x, A(:,i), k); 
        for l = 1:k
            res = res + 0.5*w(l,i)*(v(l,i)*(d_i(l)^2 + eps^2) + 1/v(l,i));
        end
    end
end

function [res] = H_epsilon(A,m,k,eps,w,x)
    setenv('distance', 'E-norm');
    res = 0;
    for i = 1:m
        d_i = distance_like(x, A(:,i), k); 
        for l = 1:k
            res = res + w(l,i)*sqrt(d_i(l)^2 + eps^2);
        end
    end
end