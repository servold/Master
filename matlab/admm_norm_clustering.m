function [X,I,t] = admm_norm_clustering(A,n,m,k,rho,max_iters,tol,X_0)
    I = zeros(1,m);
    X = X_0;
    Z = X_0;
    Y = zeros(n,k);
    ones_vec = ones(m,1);
    Phi = 0;
    weiszfeld_step_tol = 0.1;
    weiszfeld_step_max_iters = 10;
    
    for t = 1:max_iters
        prev_I = I;
        prev_Phi = Phi;
        [D,I] = clustering_distance(X, A, m, k);
        Phi = D*ones_vec;
        
        % ADMM updates for l dimension
        for l = 1:k
            CIDX_labeled_l = (I==l);
            X(:,l) = admm_weiszfeld_step(A,CIDX_labeled_l,Z(:,l),Y(:,l),rho,m,weiszfeld_step_max_iters,weiszfeld_step_tol,X(:,l));
            Z(:,l) = X(:,l) + (Y(:,l)/rho);
            Y(:,l) = Y(:,l) + rho*(X(:,l) - Z(:,l));
        end
        
        if ((sum(prev_I == I) == m) && t>1 && abs(Phi-prev_Phi)<tol)
            break;
        end
    end
    
end
