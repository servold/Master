function [f] = plot_clusters(B, I, k, x)
    plot_styles = {'b+'; 'gs'; 'ro'; 'm*';};
    f=figure();
    for l = 1:k
        CIDX_labeled_l = (I==l);
        plot(B(1,CIDX_labeled_l), B(2,CIDX_labeled_l), char(plot_styles(l)));
        hold on;
    end
    
    plot(x(1,:),x(2,:), 'kh');
    hold off;
end

