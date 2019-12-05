clc
clear

load('QS2/m10_dt0.1.mat')

s_t = [1:valMAXTIME].*valDELTIME.*2;
idx = (1/valDELTIME):length(s_t);

hFig1 = figure(1);
clf(1);
plot(s_t(idx) - 2, sum(matDGAMMADT(idx,:),2), '-ok');
grid minor
box on
axis tight

hFig2 = figure(2);
clf(2);
plot(s_t(idx) - 2, sum(matINTCIRC(idx,:),2), '-ok');
grid minor
box on
axis tight

%%
% 
% load('matlab.mat');
% s_t = [1:valMAXTIME].*valDELTIME.*2;
% idx = (1/valDELTIME):length(s_t);
% 
% hFig1 = figure(1);
% clf(1);
% plot(s_t(idx) - 2, tmp_wind(idx,3), '-ok');
% grid minor
% box on
% axis tight