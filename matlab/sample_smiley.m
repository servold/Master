function [A,n,m,k,I_true] = sample_smiley(pnts_per_eye, pnts_per_face)
    k = 3;
    n = 2;
    m = 2*pnts_per_eye + pnts_per_face;
    A = zeros(n,m);
    I_true = zeros(m,1);
    
    % smaple face
    face_center = [0.5; 0.5];
    for i = 1:pnts_per_face
        r = 0.25*sqrt(rand());
        theta = 2*pi*rand();
        A(:,i) = ([r*cos(theta); r*sin(theta)] + face_center);
        I_true(i) = 1;
    end
    
    % sample left eye
    left_eye_center = [0.2; 0.8];
    for i = (pnts_per_face+1):(pnts_per_face+pnts_per_eye)
        r = 0.05*rand();
        theta = 2*pi*rand();
        A(:,i) = ([r*cos(theta); r*sin(theta)] + left_eye_center);
        I_true(i) = 2;
    end
    
    % sample right eye
    left_eye_center = [0.8; 0.8];
    for i = (pnts_per_face+pnts_per_eye+1):m
        r = 0.05*rand();
        theta = 2*pi*rand();
        A(:,i) = ([r*cos(theta); r*sin(theta)] + left_eye_center);
        I_true(i) = 3;
    end

end

