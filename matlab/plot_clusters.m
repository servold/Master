function [] = plot_clusters(M, k, X, I, prefix)
    plot_styles = {'b+'; 'gs'; 'ro'; 'm*';};
    [~,~,iters] = size(X);
    for t = 1:iters
        if (t <=10 || mod(t,10) == 1)
            fig=figure();
            for l = 1:k
                CIDX_labeled_l = (I(:,t)==l);
                plot(M(1,CIDX_labeled_l), M(2,CIDX_labeled_l), char(plot_styles(l)));
                hold on;
            end

            plot(X(1,:,t),X(2,:,t), 'kh');
            
            ttl = ['iter-',num2str(t-1)];
            if (exist('prefix','var') == 1)
                ttl = [prefix,'-',ttl];
                saveas(fig,[ttl,'.jpg']);
            end
            title(ttl);
            hold off;
        end
    end
end

