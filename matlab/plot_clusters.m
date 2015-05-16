function [ ] = plot_clusters(B, I, k)
    plot_styles = {'b+'; 'g+'; 'r+'; 'c+'};
    figure();
    for l = 1:k
        CIDX_labeled_l = (I==l);
        plot(B(1,CIDX_labeled_l), B(2,CIDX_labeled_l), char(plot_styles(l)));
        hold on;
    end
    hold off;
end

