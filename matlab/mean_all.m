function [mean_iters,mean_time,mean_similarity,mean_psi,mean_delta_x] = mean_all(trials,max_iters,iters,time,similarity,Psi,delta_x)
    mean_iters = mean(iters,3);
    mean_time = mean(time,3);
    mean_similarity = mean(similarity,3);
    s = size(iters);
    algs = s(1);
    no_nans_psi = Psi;
    no_nans_delta_x = delta_x;
    for i = 1:algs
        for j = 1:trials
            for k = 2:max_iters
                if(isnan(no_nans_psi(i,k,j)))
                    no_nans_psi(i,k,j) = no_nans_psi(i,k-1,j);
                end
                if(isnan(no_nans_delta_x(i,k,j)))
                    no_nans_delta_x(i,k,j) = no_nans_delta_x(i,k-1,j);
                end
            end
        end
    end
    mean_psi = mean(no_nans_psi,3);
    mean_delta_x = mean(no_nans_delta_x,3);
end