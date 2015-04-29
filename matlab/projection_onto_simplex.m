function [p] = projection_onto_simplex(x)
% projection_onto_simplex computes the projection of a given vector @x of
% size 'n' onto the n-dimensional simplex
    
    if x >= 0
        if sum(x) == 1 % x is inside the simplex
            p = x;
            return
        end
    end

    n = length(x); 
    sorted_x = sort(x, 'descend'); 
    elem_sum = 0;
    delta = 0;

    % find delta:
    for i = 1:n
        elem_sum = elem_sum + sorted_x(i);
        delta = (1 - elem_sum)/i;
        
        if i == n || -delta >= sorted_x(i+1)
           break;
        end
    end

    % construct p:
    p= max(x + delta, 0); % component-wise
    
return;

