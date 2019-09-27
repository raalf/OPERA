clc
clear
% close all

A = dlmread('Alpha_15_550.txt', '', 9, 0);
A(A(:,4) > 360,:) = [];

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
load('Alpha 15 Results/relaxed_coarse.mat', 'CT_U', 'CT', 'valDELTIME', 'matDVE', 'matVLST', 'vecHUB', 'vecDVESURFACE', 'matDGAMMADT', 'matINTCIRC', 'vecTEDVE', 'tmpDVETHRUST', 'matSPANDIR')
CT_fixed = CT_U(~isnan(CT_U));
CT_fixed_s = CT(~isnan(CT));
% CT_fixed = CT_fixed(100:end);
deg_per_ts = valRPM.*(pi/30).*(180/pi).*valDELTIME;
tmp_offset = 0;
vecPOS = [tmp_offset:(length(CT_fixed)-1 + tmp_offset)]'.*deg_per_ts;

% figure(5);
% clf(5);
% plot(vecPOS_TUNNEL(100:600)./360, CT_tunnel(100:600), '-sk')
% grid minor
% box on
% axis tight


%% Plotting
hFig2 = figure(2);
clf(2);

tmp56 = [vecPOS_TUNNEL_OG CT_tunnel];

start = 75;
tmp56(:,2) = [CT_tunnel(vecPOS_TUNNEL_OG > start); CT_tunnel(vecPOS_TUNNEL_OG <= start)];
tmp56(:,1) = [vecPOS_TUNNEL(vecPOS_TUNNEL_OG > start)+start; vecPOS_TUNNEL(vecPOS_TUNNEL_OG <= start)+start];
tmp56(:,1) = rem(tmp56(:,1),360);
[a,b] = sort(tmp56(:,1));
tmp56 = [tmp56(b,:)];
tmp56(tmp56(:,2) < 3e-3,:) = [];

% scatter(tmp56(:,1), tmp56(:,2), 20, 'sk')
errorbar(tmp56(:,1), tmp56(:,2), repmat(mean(error_tunnel), size(tmp56,1), 1), 'sk')
tmp56 = smoothdata(tmp56, 'movmean',50);
hold on
plot(tmp56(:,1), tmp56(:,2),'-r','linewidth',2)
tmp = vecPOS./360;

start = 3.33;
idx = tmp >= start & tmp <= start+1;
tmp2 = find(idx);
offset = -vecPOS(tmp2(1));
plot(vecPOS(idx) + offset, CT_fixed_s(idx), '-.m');
plot(vecPOS(idx) + offset, CT_fixed(idx), '--b');

grid minor
box on
axis tight
% legend('Wind Tunnel Test Data','Moving Average of Tunnel Data','OPERA (fixed wake)','Location','SouthWest')

hold off

xt = xlim;
yt = ylim;
xlim('manual')
ylim('manual')

tmp = gca;
AR = tmp.PlotBoxAspectRatio;

%%

start = 0
loc = [start+90 0e-3; start+180 0; start+270 0e-3];

matVLST = matVLST - vecHUB;
tmp6 = acos(dot(matSPANDIR(1,:), [0 1 0], 2));
matVLST = fcnGLOBSTAR(matVLST, repmat([0 0 tmp6], size(matVLST,1), 1));
tmp7 = vecPOS(idx);
for i = 1:size(loc,1)

tmpVLST = [matVLST(:,1) matVLST(:,2)];
tmpVLST = fcnGLOBSTAR([tmpVLST tmpVLST(:,1).*0], repmat([0 0 deg2rad(-loc(i,1)+tmp7(1))], size(tmpVLST,1), 1));
tmpVLST = tmpVLST(:,1:2);
tmpVLST = ([tmpVLST(:,1).*diff(xt)./AR(1) tmpVLST(:,2).*diff(yt)./AR(2)])./3;

h1 = patch('Faces',matDVE(vecDVESURFACE == 1,:,1),'Vertices',tmpVLST + loc(i,:),'FaceColor',[153 0 0]./255,'EdgeAlpha',0.2,'FaceAlpha',0.8);
h2 = patch('Faces',matDVE(vecDVESURFACE == 2,:,1),'Vertices',tmpVLST + loc(i,:),'FaceColor',[0 0 153]./255,'EdgeAlpha',0.2,'FaceAlpha',0.8);

arrow(loc(i,:) + [30 0], loc(i,:) + [15 0])
end

legend('Wind Tunnel Test Data','Moving Average of Tunnel Data','OPERA (fixed wake, quasi-unsteady)','OPERA (fixed wake, fully unsteady)','Location','NorthWest')

xlabel('Position (Degrees)');
ylabel('C_T');


% idx = 1:500;
% CT_tunnel = CT_tunnel(idx);
% vecPOS_TUNNEL = vecPOS_TUNNEL(idx);
% 
% hFig20 = figure(20);
% clf(20);
% subplot(2,1,1)
% plot(vecPOS_TUNNEL, CT_tunnel, '-b')
% 
% offset = 907-35;
% hold on
% plot(vecPOS+offset, CT_fixed + 0.004, '-.r')
% 
% % t3 = linspace(0, valDELTIME*(length(CT_relaxed)), length(CT_relaxed)) + offset;
% % plot(t3, CT_relaxed, '-k')
% 
% grid minor
% box on
% axis tight
% hold off
% 
% subplot(2,1,2)
% yyaxis left
% Fs = 2000;            % Sampling frequency                    
% y = fft(detrend(CT_tunnel));
% f = (0:length(y)-1)*Fs/length(y);
% plot(f,abs(y)/length(y), '--b') 
% 
% hold on
% yyaxis right
% Fs = 1/valDELTIME;
% y = fft(detrend(CT_fixed));
% f = (0:length(y)-1)*Fs/length(y);
% plot(f,abs(y)/length(y),'-.r') 
% 
% % y = fft(detrend(CT_relaxed));
% % f = (0:length(y)-1)*Fs/length(y);
% % plot(f,abs(y)/length(y),'-k')
% % hold off
% 
% grid minor
% box on
% axis tight


% hFig20 = figure(20);
% clf(20);

% idx2 = vecDVESURFACE(vecTEDVE);
% yyaxis left
% plot(sum(matINTCIRC(:,idx==1),2), '-k')
% hold on
% plot(sum(matINTCIRC(:,idx==2),2), '--k')
% grid minor

% yyaxis right
% % plot(vecPOS(idx) + offset, sum(matDGAMMADT(idx,idx2==1),2), '-m')
% hold on
% % plot(vecPOS(idx) + offset, sum(matDGAMMADT(idx,idx2==2),2), '--m')
% plot(vecPOS(idx) + offset, sum(matDGAMMADT(idx,:),2), '--m')
% hold off
% % axis tight
% % grid minor
% box on

% yyaxis right
% plot(vecPOS(idx) + offset, sum(matINTCIRC(idx,idx2==1),2), '-k')
% hold on
% plot(vecPOS(idx) + offset, sum(matINTCIRC(idx,idx2==2),2), '--k')
% plot(vecPOS(idx) + offset, sum(matINTCIRC(idx,:),2), '--m')
% hold off
% axis tight
% grid minor
% box on

% yyaxis right
% tmp = vecPOS(idx) + offset;
% % plot(tmp(2:end), diff(sum(matDGAMMADT(idx,idx2==1),2)), '-k')
% hold on
% % plot(tmp(2:end), diff(sum(matDGAMMADT(idx,idx2==2),2)), '--k')
% plot(tmp(2:end), diff(sum(matDGAMMADT(idx,:),2)), '--k')
% hold off
% % axis tight
% % grid minor
% box on

% yyaxis right
% tmp = vecPOS(idx) + offset;
% plot(tmp, sum(tmpDVETHRUST(idx,idx2==1),2), '-k')
% hold on
% plot(tmp, sum(tmpDVETHRUST(idx,idx2==2),2), '--k')
% % plot(tmp(2:end), diff(sum(tmpDVETHRUST(idx,:),2)), '--k')
% hold off
% % axis tight
% % grid minor
% box on















