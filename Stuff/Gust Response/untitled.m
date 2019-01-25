clc
clear

load('steady.mat','valTIMESTEP','CL','valDELTIME')

load('Gust_CL_Time.mat')

hFig1 = figure(1);
clf(1);

plot((1:valTIMESTEP).*valDELTIME.*0.1, CL(1:valTIMESTEP),'-ok');
hold on
shift_x = 0.2;

plot(Gust_CL_Time(:,1)+shift_x, Gust_CL_Time(:,2) - Gust_CL_Time(1,2),'--^m')
hold off
box on
grid minor