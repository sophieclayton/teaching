% Make a plot of Riley's data
% Sophie Clayton, 2014
% sclayton@uw.edu

clear all

load Riley_data
figure(1);
subplot(4,1,1)
plot(YD, Ph,'-ob');ylabel('Phytoplankton');
subplot(4,1,2)
plot(YD, N,'-ob');ylabel('Nitrate');
subplot(4,1,3)
plot(YD, Tobs,'-ob');ylabel('Temperature');
subplot(4,1,4)
plot(YD, Zobs,'-ob');ylabel('Zooplankton');xlabel('Yearday')


