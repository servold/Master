function [A,n,m,k,closeness] = sample_d_spheres(pnts_per_cltr, d)
    k = 2^d;
    n = d;
    m = k*pnts_per_cltr;
    A = zeros(d,m);
    cntrs = zeros(k,d);
    closeness = ones(1,m);
    
    for l = 1:k
        cltr_cntr = 2*(de2bi(l-1,d) - 0.5*ones(1,d));
        cntrs(l,:) = cltr_cntr;
    end
    
    for l = 1:k
        cltr_cntr = cntrs(l);
        
        for i = 1:pnts_per_cltr
            inx = (l-1)*pnts_per_cltr + i;
            p = zeros(d,1);
            r = rand;
            if (d>1)
                angles = rand(d-1,1)*pi;
                angles(d-1) = angles(d-1)*2;
                for j = 1:(d-1)
                    p(j) = r*cos(angles(j));
                    r = r*sin(angles(j));
                end
            end
            p(d) = r;
            
            p = p + cltr_cntr';
            A(:,inx) = p;
            
            cntr_dis = norm(cntrs(l)-p);
            for c = 1:k
                 closeness(inx) = (closeness(inx) || (c==l || (norm(cntrs(c)-p) < cntr_dis)));
            end
        end
    end
    
end

