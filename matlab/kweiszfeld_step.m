function [x] = kweiszfeld_step(A,w,m,max_iters,tol,x_0)
    x = x_0;
    wa = zeros(m,1);
    
    for t = 1:max_iters
        for i = 1:m
            if (norm(x - A(:,i)) == 0)
                wa = zeros(m,1);
                wa(i) = 1;
                break;
            end
            wa(i) = w(i)/norm(x - A(:,i));
        end
        
        prev_x = x;
        x = A*wa/sum(wa);
        
        if (norm(prev_x - x) < tol)
            break;
        end
    end
end