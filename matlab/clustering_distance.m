% Computes two m-dimensional vectors, D(i) is the distance of a_i from its 
% closest center I(i). 
function [D,I] = clustering_distance(X, A, m, k)
   D = zeros(1,m);
   I = zeros(1,m);
   
   for i = 1:m
       a = A(:,i);
       [D(i),I(i)] = min(distance_like(X,a,k));
   end
   
return;

       