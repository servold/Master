function [X,W,d] = palm_clustering(A, v, m, k, max_iters, gamma, nu, beta, mu)

    rand_vector_ind = randperm(m);
    for l = 1:k
        X(:,l) = A(:, rand_vector_ind(l));
    end
    
    for i = 1:m
        W(:,i) = projection_onto_simplex(rand(1,k)');
    end
    
    for t = 1:max_iters
        for l = 1:k
            mean_a = (W(l,:).*v')*A';
            mean_a_sum_weight = W(l,:)*v;
            c_l = gamma(l)*max(mean_a_sum_weight, nu(l)/2);
            X(:,l) = X(:,l)*(1 - mean_a_sum_weight/c_l) + mean_a'/c_l;
        end
        
        for i = 1:m
            w = W(:,i);
            for l = 1:k
                coef = 1 - (norm(X(:,l) - A(:,i))^2)*v(i)/(beta(i)*mu(i));
                w(l) = coef * w(l);
            end
            W(:,i) = projection_onto_simplex(w);
        end
        
        d(t) = clustering_distance(X, A, v, m, k);
    end
    
return;