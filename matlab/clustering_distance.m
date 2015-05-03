function [d] = clustering_distance(X, A, v, m, k)
   d_vec = zeros(1,m);
   
   for i = 1:m
       a = A(:,i);
       d_vec(i) = min(distance_like(X,a,k));
   end
   
   d = d_vec*v;
return;

       