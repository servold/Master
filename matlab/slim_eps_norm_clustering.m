function [x,w,I,t,psi] = slim_eps_norm_clustering(A,n,m,k,max_iters,tol,x_0,w_0,eps)
    setenv('distance', 'E-norm');
    ones_vec = ones(m,1);
    alpha0 = diam(A,m);
    last_psi = 0;
    
    % x init & init clustering
    x = x_0;
    [~,CIDX] = clustering_distance(x, A, m, k);
    last_I = CIDX';
    
    % W init
    w = w_0;

    for t = 1:max_iters
        alpha = alpha0/(2^(t-1));
        
        % W update
        [v,j] = min(w*ones_vec);
        if v < 0.0000001
            disp(['cluster ', num2str(j) ,' close to zero: ', num2str(v)]);
        end
        for i = 1:m
            d = distance_like(x, A(:,i), k);
            for l = 1:k
                d(l) = sqrt(d(l)^2 + eps^2);
            end
            w(:,i) = projection_onto_simplex(w(:,i) - d/alpha);
        end
        
        % x update
        for l = 1:k
            u = zeros(m,1);
            for i = 1:m
                u(i) = w(l,i)/sqrt(norm(x(:,l) - A(:,i))^2 + eps^2);
            end
            x(:,l) = A*u/sum(u);
        end
        
        % Psi computations
        psi = 0;
        for i = 1:m
            d_i = distance_like(x, A(:,i), k);
            for l = 1:k
                psi = psi + w(l,i)*sqrt(d_i(l)^2 + eps^2);
            end
        end
        
        % clustering update
        [~,CIDX] = clustering_distance(x, A, m, k);
        I = CIDX';
        
        if ((sum(I == last_I) == m) && t>1 && (last_psi-psi)<tol)
            break;
        end
        
        last_psi = psi;
        last_I = I;
    end
end
