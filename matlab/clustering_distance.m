function [d] = clustering_distance(X, A, v, m, k)
   d = 0;
   
    for i = 1:m
       a = A(:,i);
       min_dist = inf;
       
       for l = 1:k
           curr_dist = norm(a - X(:,l));
           if (curr_dist < min_dist)
               min_dist = curr_dist;
           end
       end
       d = d + v(i)*min_dist;
   end
   
return;

       