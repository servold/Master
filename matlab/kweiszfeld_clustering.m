function [X,I,t] = kweiszfeld_clustering(A,m,k,max_iters,tol,X_0)
    I = zeros(1,m);
    X = X_0;
    ones_vec = ones(m,1);
    Phi = 0;
    weiszfeld_step_tol = 0.1;
    weiszfeld_step_max_iters = 10;
    
    for t = 1:max_iters
        prev_I = I;
        prev_Phi = Phi;
        [D,I] = clustering_distance(X, A, m, k);
        Phi = D*ones_vec;
        
        % updates for l dimension
        for l = 1:k
            CIDX_labeled_l = (I==l);
            X(:,l) = kweiszfeld_step(A,CIDX_labeled_l,m,weiszfeld_step_max_iters,weiszfeld_step_tol,X(:,l));
        end
        
        if ((sum(prev_I == I) == m) && t>1 && abs(Phi-prev_Phi)<tol)
            break;
        end
    end
    
end


