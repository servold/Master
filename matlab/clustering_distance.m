function [d,D,I] = clustering_distance(X, A, v, m, k)
   D = zeros(1,m);
   I = zeros(1,m);
   
   for i = 1:m
       a = A(:,i);
       [D(i),I(i)] = min(distance_like(X,a,k));
   end
   
   d = D*v;
return;

       