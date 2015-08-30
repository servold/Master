function [x,I,t,phi] = slim_kmeans_clustering(A,n,m,k,max_iters,tol,x_0)
    setenv('distance', 'sq-E-norm');
    
    % X init & init clustering
    x = x_0;
    [D,CIDX] = clustering_distance(x, A, m, k);
    last_I = CIDX';
    last_phi = sum(D);

    for t = 1:max_iters
        
        % X update
        for l = 1:k
            CIDX_labeled_l = (last_I==l);
            if (sum(CIDX_labeled_l) == 0)
                disp(['no points in center: ', num2str(l)]);
            end
            x(:,l) = (A*CIDX_labeled_l)/sum(CIDX_labeled_l);
        end
        
       % I update
        [D,CIDX] = clustering_distance(x, A, m, k);
        I = CIDX';
        phi = sum(D);
        
        if ((sum(I == last_I) == m) && abs(phi - last_phi)<tol)
            break;
        end
        
        last_phi = phi;
        last_I = I;
    end
end

