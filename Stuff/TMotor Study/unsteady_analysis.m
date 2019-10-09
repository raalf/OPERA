clc
clear
% close all

A = dlmread('Alpha_15_550.txt', '', 9, 0);
A(A(:,4) > 360,:) = [];
LW = 1;
%% Getting tunnel data
valDIAM = 0.4572;
valRPM = 3000;
valJ = mean(A(:,3))/((valRPM.*(pi/30)).*(valDIAM/2));

vecPOS_TUNNEL_OG = A(:,4);
vecPOS_TUNNEL = A(:,4);

idx_go = vecPOS_TUNNEL(2:end) < vecPOS_TUNNEL(1:end-1);
offset = [0 360:360:(length(find(idx_go))*360)]';
offset_idx = 1;
count = 1;
for i = 1:length(idx_go)
    vecPOS_TUNNEL(i) = vecPOS_TUNNEL(i) + offset(offset_idx);
    if idx_go(i) == true
        offset_idx = offset_idx + 1;
    end
end

vecDENSITY = (3386.39*29.23)./(287.058.*((A(:,2)-32).*(5/9) + 273.15));
CT_tunnel = A(:,7)./(vecDENSITY.*(pi.*((valDIAM/2).^2)).*(((valDIAM/2).*(valRPM.*(pi/30))).^2));
error_tunnel = (1/4)./(vecDENSITY.*(pi.*((valDIAM/2).^2)).*(((valDIAM/2).*(valRPM.*(pi/30))).^2));

%% Getting OPERA data
load('Alpha 15 Results/relaxed_fine.mat', 'CT_U', 'CT', 'valDELTIME', 'matDVE', 'matVLST', 'vecHUB', 'vecDVESURFACE', 'matDGAMMADT', 'matINTCIRC', 'vecTEDVE', 'tmpDVETHRUST', 'matSPANDIR')
CT_relaxed = CT_U(~isnan(CT_U));
CT_relaxed_s = CT(~isnan(CT));

deg_per_ts = valRPM.*(pi/30).*(180/pi).*valDELTIME;
tmp_offset = 0;
vecPOS = [tmp_offset:(length(CT_relaxed)-1 + tmp_offset)]'.*deg_per_ts;

%% Plotting
hFig2 = figure(2);
clf(2);

WH = [8 4.5];

tmp56 = [vecPOS_TUNNEL_OG CT_tunnel];

start = 0;
tmp56(:,2) = [CT_tunnel(vecPOS_TUNNEL_OG > start); CT_tunnel(vecPOS_TUNNEL_OG <= start)];
tmp56(:,1) = [vecPOS_TUNNEL(vecPOS_TUNNEL_OG > start)+start; vecPOS_TUNNEL(vecPOS_TUNNEL_OG <= start)+start];
tmp56(:,1) = rem(tmp56(:,1),360);
[a,b] = sort(tmp56(:,1));
tmp56 = [tmp56(b,:)];
tmp56(tmp56(:,2) < 3e-3,:) = [];

scatter(tmp56(:,1), tmp56(:,2), 5, 'sk')


% hFig40 = figure(40);
% clf(40);
% Fs = 2000;                      % samples per second
% dt = 1/Fs;                     % seconds per sample
% StopTime = 0.02;                  % seconds
% t = (0:dt:StopTime-dt)';
% N = size(t,1);
% x = linspace(0,360,length(t));
% [tmp57,idx] = uniquetol(tmp56(:,1), 1e-2);
% tmp57 = [tmp57, tmp56(idx,2)];
% y = interp1(tmp57(:,1), tmp57(:,2), x', 'linear','extrap');
% y = fft(detrend(y));
% dF = Fs/N;
% f = -Fs/2:dF:Fs/2-dF;
% plot(f,abs(y)/N, '-k')

% hFig40 = figure(40);
% clf(40);
% Fs = 2000;                      % samples per second
% dt = 1/Fs;                     % seconds per sample
% StopTime = length(CT_tunnel)*dt;                  % seconds
% t = (0:dt:StopTime-dt)';
% N = size(t,1);
% % x = linspace(0,360,length(t));
% % [tmp57,idx] = uniquetol(tmp56(:,1), 1e-2);
% % tmp57 = [tmp57, tmp56(idx,2)];
% % y = interp1(tmp57(:,1), tmp57(:,2), x', 'linear','extrap');
% y = fft(detrend(CT_tunnel));
% dF = Fs/N;
% f = -Fs/2:dF:Fs/2-dF;
% plot(f,abs(y)/N, '-k')

%%
figure(2);
tmp56 = smoothdata(tmp56, 'movmean',50);
hold on
plot(tmp56(:,1), tmp56(:,2),'-r','linewidth',LW)
tmp = vecPOS./360;

start = 2.0;
idx = tmp >= start & tmp < start+1;
idx2 = tmp >= start+1 & tmp < start+2;
idx2 = idx;
tmp2 = find(idx);
offset = -vecPOS(tmp2(1));
plot(vecPOS(idx) + offset, (CT_relaxed(idx) + CT_relaxed(idx2))./2, '--b','linewidth',LW);

% figure(40);
% dt = valDELTIME;                     % seconds per sample
% Fs = 1/dt;                      % samples per second
% StopTime = 1/(valRPM/60);                  % seconds
% t = (0:dt:StopTime-dt)';
% N = size(t,1);
% y = fft(detrend((CT_relaxed(idx) + CT_relaxed(idx2))./2));
% dF = Fs/N;
% f = -Fs/2:dF:Fs/2-dF;
% hold on
% plot(f,abs(y)/N, '--b')
% hold off

%%
% hfig5 = figure(5)
% clf(5)
% 
% plot(CT_relaxed, '-k');
% hold on
% plot(CT_relaxed_s, '--b');
% hold off
% grid minor
% box on
% axis tight
% 
% CT_relaxed = CT_U(~isnan(CT_U));
% CT_relaxed_s = CT(~isnan(CT));


%%
hFig20 = figure(20);
clf(20);

WH = [4.5 5];
idx3 = vecDVESURFACE(vecTEDVE);
tmp5 = vecPOS(idx) + offset;
plot(tmp5, (sum(tmpDVETHRUST(idx,idx3==1),2) + sum(tmpDVETHRUST(idx2,idx3==1),2))./2, '-.r','linewidth',LW)
hold on
plot(tmp5, (sum(tmpDVETHRUST(idx,idx3==2),2) + sum(tmpDVETHRUST(idx2,idx3==2),2))./2, '--b','linewidth',LW)
plot(tmp5, (sum(tmpDVETHRUST(idx,:),2) + sum(tmpDVETHRUST(idx2,:),2))./2, '-k','linewidth',LW)
hold off
axis tight
grid minor
box on

xlabel('Position (Degrees)');
ylabel('Thrust (N)');
legend('Blade 1','Blade 2', 'Combined','Location','NorthWest','FontSize',10)

xtf = xlim;
ytf = ylim;

%%
figure(2);
load('Alpha 15 Results/fixed_fine.mat', 'CT_U', 'CT', 'valDELTIME', 'matDVE', 'matVLST', 'vecHUB', 'vecDVESURFACE', 'matDGAMMADT', 'matINTCIRC', 'vecTEDVE', 'tmpDVETHRUST', 'matSPANDIR')
CT_fixed = CT_U(~isnan(CT_U));
CT_fixed_s = CT(~isnan(CT));
% CT_fixed = CT_fixed(100:end);
deg_per_ts = valRPM.*(pi/30).*(180/pi).*valDELTIME;
tmp_offset = 0;
vecPOS = [tmp_offset:(length(CT_fixed)-1 + tmp_offset)]'.*deg_per_ts;

grid minor
box on
axis tight

hold off

tmp = vecPOS./360;
start = 3.0;
idx = tmp >= start & tmp < start+1;
idx2 = tmp >= start+1 & tmp <= start+2;
tmp2 = find(idx);
offset = -vecPOS(tmp2(1));
hold on
plot(vecPOS(idx) + offset, (CT_fixed(idx) + CT_fixed(idx2))./2, '-.m','linewidth',LW);
hold off

legend('RU Test Data','Moving Average of Test Data','DDE Method (relaxed wake)','DDE Method (fixed wake)','Location','South','FontSize',8)
xlabel('Position (Degrees)');
ylabel('C_T');
WH = [4.5*2 5];
% fcnFIG2LATEX(hFig2, 'tmotor_time.pdf', WH)

% figure(40);
% dt = valDELTIME;                     % seconds per sample
% Fs = 1/dt;                      % samples per second
% StopTime = 1/(valRPM/60);                  % seconds
% t = (0:dt:StopTime-dt)';
% N = size(t,1);
% y = fft(detrend((CT_fixed(idx) + CT_fixed(idx2))./2));
% dF = Fs/N;
% f = -Fs/2:dF:Fs/2-dF;
% hold on
% plot(f,abs(y)/N, '-.m')
% hold off

%%
hFig21 = figure(21);
clf(21);

WH = [4.5 5];

idx3 = vecDVESURFACE(vecTEDVE);
tmp5 = vecPOS(idx) + offset;
plot(tmp5, (sum(tmpDVETHRUST(idx,idx3==1),2) + sum(tmpDVETHRUST(idx2,idx3==1),2))./2, '-.r','linewidth',LW)
hold on
plot(tmp5, (sum(tmpDVETHRUST(idx,idx3==2),2) + sum(tmpDVETHRUST(idx2,idx3==2),2))./2, '--b','linewidth',LW)
plot(tmp5, (sum(tmpDVETHRUST(idx,:),2) + sum(tmpDVETHRUST(idx2,:),2))./2, '-k','linewidth',LW)
hold off
axis tight
grid minor
box on

xlabel('Position (Degrees)');
ylabel('Thrust (N)');
legend('Blade 1','Blade 2', 'Combined','Location','NorthWest','FontSize',10)

% fcnFIG2LATEX(hFig21, 'blade_thrust_fixed.pdf', WH)

xtr = xlim;
ytr = ylim;

figure(20);
xlim([min([xtf(1) xtr(1)]) max([xtf(2) xtr(2)])]);
ylim([min([ytf(1) ytr(1)]) max([ytf(2) ytr(2)])]);
% fcnFIG2LATEX(hFig20, 'blade_thrust_relaxed.pdf', WH)
