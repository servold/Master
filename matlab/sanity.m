% Sanity check for PALM clustering and comparison with KMEANS clustering.
% Generates 4 clusters in R^2, with centers in @cluster_center, and the
% reflected points over the axes. Each cluster has @points_in_cluster.
% PALM algorithm runs at most @max_iters iterations. The results are 
% returned in C_palm and I_palm.
function [A,C_palm,I_palm,C_kmeans,I_kmeans] = sanity(points_in_cluster, cluster_center, max_iters, tol)
    k = 4;
    n = 2;
    m = k*points_in_cluster;
    A = zeros(n,m);
    cluster_center = cluster_center - [1,1];
    
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
    X_0 = [p1'+[0.1,0.1]',p1'+[0.1,-0.1]',p1'+[-0.1,0.1]',p3'+[-0.1,-0.1]'];
    disp('finished precomputing');
    
    %PALM clustering
    tic;
    [X_pc, I_pc, palm_iters] = palm_clustering(A,n,m,k,max_iters,tol,X_0);
    C_palm = X_pc(:,:,palm_iters+1);
    I_palm = I_pc(:,palm_iters+1);
    disp(['PALM clustering time: ', num2str(toc), ', # iterations: ', num2str(palm_iters)]);
    
    
    %KMEANS clustering
    tic;
    [X_kmc, I_kmc, kmeans_iters] = kmeans_clustering(A,n,m,k,max_iters,tol,X_0);
    C_kmeans = X_kmc(:,:,kmeans_iters+1);
    I_kmeans = I_kmc(:,kmeans_iters+1);
    disp(['KMEANS clustering time: ', num2str(toc), ', # iterations: ', num2str(kmeans_iters)]);
    
    plot_clusters(A, k, X_pc, I_pc);
    plot_clusters(A, k, X_kmc, I_kmc);
    
end

