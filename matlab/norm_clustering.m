function [x,I,t] = norm_clustering(A,n,m,k,max_iters,tol,x_0,prefix)
    w = zeros(k,m);
    Psi = zeros(1,max_iters);
    ones_vec = ones(m,1);
    alpha = 10;
    
    % x init & init clustering
    x = x_0;
    [~,CIDX] = clustering_distance(x, A, m, k);
    I = CIDX';
    
    % plot x0
    fig = plot_clusters(A,I,k,x);
%    saveas(fig,[prefix,'-alpha-',num2str(alpha),'-iter-',num2str(0),'.jpg']);
    
    % W init
    for i = 1:m
        w(:,i) = rand(k,1);
        w(:,i) = w(:,i)/sum(w(:,i));
    end

    for t = 1:max_iters
        alpha = 100/t;
        
        % w update
        [v,j] = min(w*ones_vec);
        if v < 0.0000001
            disp(['cluster ', num2str(j) ,' close to zero: ', num2str(v)]);
        end
        for i = 1:m
            d = distance_like(x, A(:,i), k);
            w(:,i) = projection_onto_simplex(w(:,i) - d/alpha);
        end
        
        % x update
        for l = 1:k
            x_l = x(:,l);
            w_l = w(l,:);
            u = zeros(m,1);
            
            for i = 1:m
                u(i) = w_l(i)/norm(x_l - A(:,i));
            end
            
            x(:,l) = A*u/sum(u);
        end
        
        % Psi computations
        for i = 1:m
            Psi(t) = Psi(t) + distance_like(x, A(:,i), k)'*w(:,i);
        end
        
        % clustering update
        [~,CIDX] = clustering_distance(x, A, m, k);
        I_new = CIDX';
        
        if (t <=10 || mod(t,10) == 1)
            fig = plot_clusters(A,I_new,k,x);
%            saveas(fig,[prefix,'-alpha-',num2str(alpha),'-iter-',num2str(t),'.jpg']);
        end
        
        if ((sum(I == I_new) == m) && t>1 && (Psi(t-1)-Psi(t))<tol)
            break;
        end
        I = I_new;
    end
    
end
