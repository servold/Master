function [A,n,m,k,I_true] = sample_ellipse(num_on_points)
    k = 2;
    n = 2;
    m = num_on_points;
    A = zeros(n,m);
    I_true = zeros(m,1);
    
    a = 0.2;
    b = 0.3;
    left_center = [0.5;0.6];
    right_center = [0.5;0.3];
    for i = 1:m
        theta = 2*pi*rand();
        if (theta > pi/2 & theta < 3*pi/2)
            I_true(i) = 1;
            A(:,i) = ([a*cos(theta); b*sin(theta)] + left_center);
        else
            I_true(i) = 2;
            A(:,i) = ([a*cos(theta); b*sin(theta)] + right_center);
        end
    end

end


