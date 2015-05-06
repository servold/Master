function [renamed_clusters] = rename_clusters(clusters, naming_perm, k)
    renamed_clusters = zeros(1,k);
    
    for l = 1:k
       renamed_clusters(clusters == l) = naming_perm(l);
    end
end

