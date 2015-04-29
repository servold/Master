function [r] = gradient_steps_ratio(n, m, k, max_iters)
    r = zeros(1,max_iters*m);
    A = zeros(n,m);
    X = zeros(n,k);
    Y = zeros(n,k);
    dX = zeros(1,k);
    dY = zeros(1,k);
    j = 1;

    for t = 1:max_iters
        for i = 1:m
            A(:,i) = randn(n,1);
        end
        
        for l = 1:k
            X(:,l) = randn(n,1);
            Y(:,l) = randn(n,1);
        end
        norm_X_Y = norm(X(:)-Y(:));
        
        for i = 1:m
            for l = 1:k
                dX(l) = mpower(norm(X(:,l)-A(:,i)),2);
                dY(l) = mpower(norm(Y(:,l)-A(:,i)),2);
            end
            r(j) = (norm_X_Y/norm(dX-dY));
            j = j+1;
        end
    end
end