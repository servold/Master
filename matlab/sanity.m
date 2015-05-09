% Sanity check for PALM clustering and comparison with KMEANS clustering.
% Generates 4 clusters in R^2, with centers in @cluster_center, and the
% reflected points over the axes. Each cluster has @points_in_cluster.
% PALM algorithm runs for @iters iterations. The results are returned in
% C_palm and I_palm, S_palm is the starting point of the algorithm.
function [A,S_palm,C_palm,I_palm,Psi,C_kmeans,I_kmeans] = sanity(points_in_cluster, cluster_center, iters)
    k = 4;
    n = 2;
    m = k*points_in_cluster;
    A = zeros(n,m);
    cluster_center = cluster_center - [1,1];
    plot_styles = {'b+'; 'g+'; 'r+'; 'c+'};
    
    for i = 1:points_in_cluster
        p1 = (cluster_center + [rand()*2, rand()*2]).*[+1,+1];
        p2 = (cluster_center + [rand()*2, rand()*2]).*[-1,+1];
        p3 = (cluster_center + [rand()*2, rand()*2]).*[-1,-1];
        p4 = (cluster_center + [rand()*2, rand()*2]).*[+1,-1];
        A(:, (i-1)*k + 1) = p1';
        A(:, (i-1)*k + 2) = p2';
        A(:, (i-1)*k + 3) = p3';
        A(:, (i-1)*k + 4) = p4';
    end
    
    %PALM clustering
    [X,W,Psi] = palm_clustering(A,n,m,k,iters);
    S_palm = X(:,:,1);
    C_palm = X(:,:,iters+1);
    [d,D,I_palm] = clustering_distance(C_palm, A, ones(m,1), m, k);
    
    figure();
    for i = 1:m
        plot(A(1,i), A(2,i), char(plot_styles(I_palm(i))));
        hold on;
    end
    title('PALM clustering');
    hold off;
    
    %KMEANS clustering
    [I_kmeans, C_kmeans] = kmeans(A', k);
    figure();
    for i = 1:m
        plot(A(1,i), A(2,i), char(plot_styles(I_kmeans(i))));
        hold on;
    end
    title('KMEANS clustering');
    hold off;
end

