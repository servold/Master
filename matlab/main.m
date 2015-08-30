function [R] = main(trials, A, n, m, k, max_iters, tol, eps)
    R = zeros(4,3,trials);
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
        [~,~,t,psi] = slim_palm_clustering(A,n,m,k,max_iters,tol,x0,w0);
        y = toc;
        R(3,1,j) = t;
        R(3,2,j) = y;
        R(3,3,j) = psi;
        
        tic;
        [~,~,t,psi] = slim_eps_norm_clustering(A,n,m,k,max_iters,tol,x0,w0,eps);
        y = toc;
        R(4,1,j) = t;
        R(4,2,j) = y;
        R(4,3,j) = psi;
    end
end