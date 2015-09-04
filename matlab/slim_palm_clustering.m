function [x,I,t,psi] = slim_palm_clustering(A,n,m,k,max_iters,tol,x_0,w_0)
    setenv('distance', 'sq-E-norm');
    ones_vec = ones(m,1);
    alpha0 = diam(A,m);
    last_psi = 0;
    
    % X init & init clustering
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
            w(:,i) = projection_onto_simplex(w(:,i) - d/alpha);
        end
        
        % X update
        for l = 1:k
            x(:,l) = (A*w(l,:)')/sum(w(l,:));
        end
        
        psi = 0;
        % Psi computations
        for i = 1:m
            psi = psi + distance_like(x, A(:,i), k)'*w(:,i);
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

