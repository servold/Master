function [iters,time,similarity,Psi,delta_x] = main(trials, A, n, m, k, max_iters, tol, eps, I_true)
%     alpha_update_functions = {@(a,t)a/(2^(t-1)), @(a,t)a/t,
%     @(a,t)a/(t^2), @(a,t)a, @(a,t)10*a, @(a,t)0.1*a};
    alpha_update_functions = {@(a,t)a/(2^(t-1)), @(a,t)a/(t^2), @(a,t)a, @(a,t)0.1*a};
    iters = zeros(2+3*length(alpha_update_functions),1,trials);
    time = zeros(2+3*length(alpha_update_functions),1,trials);
    similarity = zeros(2+3*length(alpha_update_functions),1,trials);
    Psi = zeros(2+3*length(alpha_update_functions),max_iters,trials);
    Psi(:,:,:)=nan;
    delta_x = zeros(2+3*length(alpha_update_functions),max_iters-1,trials);
    delta_x(:,:,:)=nan;
    for j = 1:trials
        disp(['trial ', num2str(j)]);
        [w0,x0] = rand_init(A,n,k,m);
        
        tic;
        [x,I,t,psi] = kmeans_clustering(A,n,m,k,max_iters,tol,x0);
        y = toc;
        iters(1,1,j) = t;
        time(1,1,j) = y;
        similarity(1,1,j) = compare_clusters(I_true, I(:,t+1), k);
        Psi(1,1:(t+1),j) = psi(1:t+1);
        for i = 2:t
            delta_x(1,i,j) = norm(x(:,:,i)-x(:,:,i-1))/norm(x(:,:,i-1));
        end
        
        tic;
        x_pp = kmeans_pp_init(A,n,m,k);
        [x,I,t,psi] = kmeans_clustering(A,n,m,k,max_iters,tol,x_pp);
        y = toc;
        iters(2,1,j) = t;
        time(2,1,j) = y;
        similarity(2,1,j) = compare_clusters(I_true, I(:,t+1), k);
        Psi(2,1:(t+1),j) = psi(1:t+1);
        for i = 2:t
            delta_x(2,i,j) = norm(x(:,:,i)-x(:,:,i-1))/norm(x(:,:,i-1));
        end
        
        for f = 1:length(alpha_update_functions)
            idx = 2+f;
            tic;
            [x,w,~,I_W,t,psi] = palm_clustering(A,n,m,k,max_iters,tol,x0,w0,alpha_update_functions{f});
            y = toc;
            iters(idx,1,j) = t;
            time(idx,1,j) = y;
            similarity(idx,1,j) = compare_clusters(I_true, I_W(:,t+1), k);
            Psi(idx,1:length(psi),j) = psi;
            for i = 2:t
                delta_x(idx,i,j) = norm(x(:,:,i)-x(:,:,i-1))/norm(x(:,:,i-1));
            end
        end
        
        for f = 1:length(alpha_update_functions)
            idx = 2+length(alpha_update_functions)+f;
            tic;
            [x,w,~,I_W,t,psi] = palm_clustering(A,n,m,k,max_iters,tol,x_pp,w0,alpha_update_functions{f});
            y = toc;
            iters(idx,1,j) = t;
            time(idx,1,j) = y;
            similarity(idx,1,j) = compare_clusters(I_true, I_W(:,t+1), k);
            Psi(idx,1:length(psi),j) = psi;
            for i = 2:t
                delta_x(idx,i,j) = norm(x(:,:,i)-x(:,:,i-1))/norm(x(:,:,i-1));
            end
        end

%         for f = 1:length(alpha_update_functions)
%             idx = f;
%             tic;
%             [x,w,~,I_W,t,psi] = eps_norm_clustering(A,n,m,k,max_iters,tol,x0,w0,eps,alpha_update_functions{f});
%             y = toc;
%             iters(idx,1,j) = t;
%             time(idx,1,j) = y;
%             similarity(idx,1,j) = compare_clusters(I_true, I_W(:,t+1), k);
%             Psi(idx,1:length(psi),j) = psi;
%             for i = 2:t
%                 delta_x(idx,i,j) = norm(x(:,:,i)-x(:,:,i-1))/norm(x(:,:,i-1));
%             end
%         end
        
        for f = 1:length(alpha_update_functions)
            idx = 2+2*length(alpha_update_functions)+f;
            tic;
            [x,w,~,~,~,I_W,t,~,psi,~,~] = eps_norm_clustering3(A,n,m,k,max_iters,tol,x0,w0,eps,alpha_update_functions{f});
            y = toc;
            iters(idx,1,j) = t;
            time(idx,1,j) = y;
            similarity(idx,1,j) = compare_clusters(I_true, I_W(:,t+1), k);
            Psi(idx,1:length(psi),j) = psi;
            for i = 2:t
                delta_x(idx,i,j) = norm(x(:,:,i)-x(:,:,i-1))/norm(x(:,:,i-1));
            end
        end

    end
end