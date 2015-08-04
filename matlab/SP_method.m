function [w0,x0] = SP_method(A,n,k,m)
    w0 = zeros(k,m);
    for i = 1:m
        w0(:,i) = rand(k,1);
        w0(:,i) = w0(:,i)/sum(w0(:,i));
    end
    
    x0 = zeros(n,k);
    for l = 1:k
        w0_l = w0(l,:);
        
        pl = 1;
        min_hl = H_l(w0_l, A(:,1), A,m);
        for i = 2:m
            if H_l(w0_l, A(:,i), A,m) < min_hl
                pl = i;
                min_hl = H_l(w0_l, A(:,i), A,m);
            end
        end
        
        xl = A(:,pl);
        Rl = 0;
        for i = 1:m
            if i~= pl
                Rl = Rl + w0_l(i)*(xl - A(:,i))/norm(xl - A(:,i));
            end
        end
        
        if norm(Rl) ~= 0
            Ll = 0;
            for i = 1:m
                if i~= pl
                    Ll = Ll + w0_l(i)/norm(xl - A(:,i));
                end
            end
            
            tl = (norm(Rl) - w0_l(pl))/Ll;
            dl = -Rl/norm(Rl);
            x0(:,l) = A(:,pl) + tl*dl;
        else
            x0(:,l)=A(:,pl);
        end
    end
    
end