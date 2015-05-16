function [ ] = plot_clusters(B, I, m)
    plot_styles = {'b+'; 'g+'; 'r+'; 'c+'};
    figure();
    for i = 1:m
        plot(B(1,i), B(2,i), char(plot_styles(I(i))));
        hold on;
    end
    hold off;
end

