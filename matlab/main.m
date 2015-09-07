function [R] = main(trials, A, n, m, k, max_iters, tol, eps)
    R = zeros(4,4,trials);
    for j = 1:trials
        [w0,x0] = rand_init(A,n,k,m);
        
        tic;
        [x,I,t,phi] = slim_kmeans_clustering(A,n,m,k,max_iters,tol,x0);
        y = toc;
        R(1,1,j) = t;
        R(1,2,j) = y;
        R(1,3,j) = phi;
        
        tic;
        [I,X,t] = kmeans(A,k);
        y = toc;
        R(2,1,j) = t;
        R(2,2,j) = y;
        [D,I] = clustering_distance(X, A, m, k);
        R(2,3,j) = sum(D);
        
        tic;
        [x,w,~,t,psi] = slim_palm_clustering(A,n,m,k,max_iters,tol,x0,w0);
        y = toc;
        R(3,1,j) = t;
        R(3,2,j) = y;
        R(3,3,j) = psi;
        
        setenv('distance', 'E-norm');
        % Psi computations
        comp_psi = 0;
        for i = 1:m
            d_i = distance_like(x, A(:,i), k);
            for l = 1:k
                comp_psi = comp_psi + w(l,i)*sqrt(d_i(l)^2 + eps^2);
            end
        end
        R(3,4,j) = comp_psi;
        
        tic;
        [x,w,~,t,psi] = slim_eps_norm_clustering(A,n,m,k,max_iters,tol,x0,w0,eps);
        y = toc;
        R(4,1,j) = t;
        R(4,2,j) = y;
        R(4,3,j) = psi;
        
        setenv('distance', 'sq-E-norm');
        comp_psi = 0;
        % Psi computations
        for i = 1:m
            comp_psi = comp_psi + distance_like(x, A(:,i), k)'*w(:,i);
        end
        R(4,4,j) = comp_psi;
    end
end