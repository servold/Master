function [iters,time,Psi,delta_x] = main(trials, A, n, m, k, max_iters, tol, eps)
    iters = zeros(4,1,trials);
    time = zeros(4,1,trials);
    Psi = zeros(4,max_iters,trials);
    delta_x = zeros(4,max_iters-1,trials);
    for j = 1:trials
        [w0,x0] = rand_init(A,n,k,m);
        
        tic;
        [x,I,t,psi] = kmeans_clustering(A,n,m,k,max_iters,tol,x0);
        y = toc;
        iters(1,1,j) = t;
        time(1,1,j) = y;
        Psi(1,:,j) = psi(2:t+1);
        for i = 2:t
            delta_x(1,i,j) = norm(x(:,:,i)-x(:,:,i-1))/norm(x(:,:,i-1));
        end
        
        tic;
        x_pp = kmeans_pp_init(A,n,m,k);
        [x,I,t,psi] = kmeans_clustering(A,n,m,k,max_iters,tol,x_pp);
        y = toc;
        iters(2,1,j) = t;
        time(2,1,j) = y;
        Psi(2,:,j) = psi(2:t+1);
        for i = 2:t
            delta_x(2,i,j) = norm(x(:,:,i)-x(:,:,i-1))/norm(x(:,:,i-1));
        end
        
        tic;
        [x,w,~,t,psi] = palm_clustering(A,n,m,k,max_iters,tol,x0,w0);
        y = toc;
        iters(3,1,j) = t;
        time(3,1,j) = y;
        Psi(3,:,j) = psi;
        for i = 2:t
            delta_x(3,i,j) = norm(x(:,:,i)-x(:,:,i-1))/norm(x(:,:,i-1));
        end
        
%         setenv('distance', 'E-norm');
%         % Psi computations
%         comp_psi = 0;
%         for i = 1:m
%             d_i = distance_like(x, A(:,i), k);
%             for l = 1:k
%                 comp_psi = comp_psi + w(l,i)*sqrt(d_i(l)^2 + eps^2);
%             end
%         end
%         R(3,4,j) = comp_psi;
        
        tic;
        [x,w,~,t,psi] = eps_norm_clustering(A,n,m,k,max_iters,tol,x0,w0,eps);
        y = toc;
        iters(4,1,j) = t;
        time(4,1,j) = y;
        Psi(4,:,j) = psi;
        for i = 2:t
            delta_x(4,i,j) = norm(x(:,:,i)-x(:,:,i-1))/norm(x(:,:,i-1));
        end
        
%         setenv('distance', 'sq-E-norm');
%         comp_psi = 0;
%         % Psi computations
%         for i = 1:m
%             comp_psi = comp_psi + distance_like(x, A(:,i), k)'*w(:,i);
%         end
%         R(4,4,j) = comp_psi;
    end
end