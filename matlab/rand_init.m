function [w0,x0] = SP_method(A,n,k,m)
    w0 = zeros(k,m);
    for i = 1:m
        w0(:,i) = rand(k,1);
        w0(:,i) = w0(:,i)/sum(w0(:,i));
    end
    
    rand_idx = randperm(m);
    x0 = A(:,rand_idx(1:k));
end