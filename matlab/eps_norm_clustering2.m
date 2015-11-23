function [X,W,V,I,I_W,t,Psi] = eps_norm_clustering2(A,n,m,k,max_iters,tol,x_0,w_0,eps,alpha_update_f)
    setenv('distance', 'E-norm');
    X = zeros(n,k,(max_iters+1));
    W = zeros(k,m,(max_iters+1));
    V = zeros(k,m,(max_iters+1));
    I = zeros(m,(max_iters+1));
    I_W = zeros(m,(max_iters+1));
    Psi = zeros(1,max_iters+1);
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
    for j = 1:m
        [~,b] = max(w_0(:,j));
        I_W(j,1) = b;
    end
    
    % V init
    V(:,:,1) = mean([1/kapa, 1/eps])*ones(k,m);
    
    % Psi computations
    for i = 1:m
        d_i = distance_like(X(:,:,1), A(:,i), k);
        w_i = W(:,i,1);
        v_i = V(:,i,1);
        for l = 1:k
            Psi(1) = Psi(1) + 0.5*w_i(l)*(v_i(l)*(d_i(l)^2 + eps^2) + w_i(l)/v_i(l));
        end
    end

    for t = 1:max_iters
        alpha = max(1e-8,alpha_update_f(alpha0,t));
        
        % W&V update
        [v,j] = min(W(:,:,t)*ones_vec);
        if v < 0.0000001
            disp(['cluster ', num2str(j) ,' close to zero: ', num2str(v)]);
        end
        for i = 1:m
            d_i = distance_like(X(:,:,t), A(:,i), k);
            for l = 1:k
                V(l,i,t+1) = 1/sqrt(d_i(l)^2 + eps^2);
                d_i(l) = 0.5*(V(l,i,t)*(d_i(l)^2 + eps^2) + 1/V(l,i,t));
            end
            W(:,i,t+1) = projection_onto_simplex(W(:,i,t) - d_i/alpha);
        end
        
        % X update
        for l = 1:k
            WV_l_row = W(l,:,t+1).*V(l,:,t+1);
            X(:,l,t+1) = (A*WV_l_row')/sum(WV_l_row);
        end
        
        % Psi computations
        for i = 1:m
            d_i = distance_like(X(:,:,t+1), A(:,i), k);
            w_i = W(:,i,t+1);
            v_i = V(:,i,t+1);
            for l = 1:k
                Psi(t+1) = Psi(t+1) + 0.5*w_i(l)*(v_i(l)*(d_i(l)^2 + eps^2) + w_i(l)/v_i(l));
            end
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
    V = V(:,:,1:t+1);
    I = I(:,1:t+1);
    I_W = I_W(:,1:t+1);
    Psi = Psi(:,1:t+1);
end

