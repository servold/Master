function [ ] = plot_results(mean_psi)
    fig=figure();
    plot(mean_psi([1,2,3,9],:)')
    legend('kmeans','kmeans++','kpalm, \alpha(t)=diam(A)/2^{t-1}','kpalm++, \alpha(t)=diam(A)/2^{t-1}')
    saveas(fig,'1.fig');
    fig=figure();
    plot(mean_psi([1,2,4,10],:)')
    legend('kmeans','kmeans++','kpalm, \alpha(t)=diam(A)/t','kpalm++, \alpha(t)=diam(A)/t')
    saveas(fig,'2.fig');
    fig=figure();
    plot(mean_psi([1,2,5,11],:)')
    legend('kmeans','kmeans++','kpalm, \alpha(t)=diam(A)/t^2','kpalm++, \alpha(t)=diam(A)/t^2')
    saveas(fig,'3.fig');
    fig=figure();
    plot(mean_psi([1,2,6,12],:)')
    legend('kmeans','kmeans++','kpalm, \alpha(t)=diam(A)','kpalm++, \alpha(t)=diam(A)')
    saveas(fig,'4.fig');
    fig=figure();
    plot(mean_psi([1,2,7,13],:)')
    legend('kmeans','kmeans++','kpalm, \alpha(t)=10diam(A)','kpalm++, \alpha(t)=diam(A)')
    saveas(fig,'5.fig');
    fig=figure();
    plot(mean_psi([1,2,8,14],:)')
    legend('kmeans','kmeans++','kpalm, \alpha(t)=0.1diam(A)','kpalm++, \alpha(t)=0.1diam(A)')
    saveas(fig,'6.fig');
    fig=figure();
    plot(mean_psi(3:8,:)')
    legend('\alpha(t)=diam(A)/2^{t-1}','\alpha(t)=diam(A)/t','\alpha(t)=diam(A)/t^2','\alpha(t)=diam(A)','\alpha(t)=10diam(A)','\alpha(t)=0.1diam(A)')
    saveas(fig,'7.fig');
    fig=figure();
    plot(mean_psi(9:14,:)')
    legend('\alpha(t)=diam(A)/2^{t-1}','\alpha(t)=diam(A)/t','\alpha(t)=diam(A)/t^2','\alpha(t)=diam(A)','\alpha(t)=10diam(A)','\alpha(t)=0.1diam(A)')
    saveas(fig,'8.fig');
end

