function [d] = clustering_distance(X, A, v, m, k)
   d = 0;
   
    for i = 1:m
       a = A(:,i);
       d = d + v(i)*min(distance_like(X,a));
   end
   
return;

       