function [u] = admm_weiszfeld_step(A,w,z,y,rho,m,max_iters,tol,u_0)
    u = u_0;
    wa = zeros(m,1);
    
    for t = 1:max_iters
        for i = 1:m
            wa(i) = w(i)/norm(u - A(:,i));
        end
        
        prev_u = u;
        u = (rho*z - y + A*wa)/(rho + sum(wa));
        
        if (norm(prev_u - u) < tol)
            break;
        end
    end
end

