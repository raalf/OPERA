clc
clear

load('G:\GIT\opera\Stuff\TMotor Study\Tunnel Testing\2019-10-31\31-Oct-2019 12.41.35_Scorpion_ASI_T-Motor 18in_RPM3000_Alpha30_21.5736.mat') % 0.3004
% load('G:\GIT\opera\Stuff\TMotor Study\Tunnel Testing\2019-11-06\06-Nov-2019 12.34.15_Scorpion_ASI_T-Motor 18in_RPM3000_Alpha0_9.667.mat') % 0.1346

close all

%% Getting OPERA data
load('Alpha 0 Results/TMotor_Relaxed_J0.3_0.0005.mat', 'CT_U', 'CT', 'valDELTIME', 'matDVE', 'matVLST', 'vecHUB', 'vecDVESURFACE', 'matDGAMMADT', 'matINTCIRC', 'vecTEDVE', 'valMAXTIME', 'matSPANDIR', 'valRPM', 'vecDVETHRUST')
% load('Alpha 0 Results/TMotor_Relaxed_J0.1346.mat', 'CT_U', 'CT', 'valDELTIME', 'matDVE', 'matVLST', 'vecHUB', 'vecDVESURFACE', 'matDGAMMADT', 'matINTCIRC', 'vecTEDVE', 'vecDVETHRUST', 'matSPANDIR', 'valRPM')

CT_relaxed = CT_U(~isnan(CT_U));
% CT_relaxed = CT(~isnan(CT));

deg_per_ts = valRPM.*(pi/30).*(180/pi).*valDELTIME;
tmp_offset = 0;
vecPOS_R = [tmp_offset:(length(CT_relaxed)-1 + tmp_offset)]'.*deg_per_ts;

%% Plotting
hFig2 = figure(2);
clf(2);

binAng = linspace(0, 360, 30);
for i = 1:length(binAng)
    rng = 0.25;
    idx = vecPOS_TUNNEL_OG >= binAng(i) - rng & vecPOS_TUNNEL_OG <= binAng(i) + rng;
    
    binAvg(i) = mean(CT_tunnel(idx));
    binMax(i) = max(CT_tunnel(idx));
    binMin(i) = min(CT_tunnel(idx));
end
errorbar(binAng, binAvg, binMin - binAvg, binMax - binAvg, ':sk');

%%

tmp = ((vecPOS_R + 90)./360);
start = 2;
idx = tmp >= start & tmp < start+1;
idx2 = tmp >= start+1 & tmp < start+2;
% idx2 = idx;
tmp2 = find(idx);
offset = -vecPOS_R(tmp2(1));
hold on
plot(vecPOS_R(idx) + offset, (CT_relaxed(idx) + CT_relaxed(idx2))./2, '--b');
hold off
grid minor
box on
axis tight

% tmp = ((vecPOS_F + 90)./360);
% start = 2;
% idx = tmp >= start & tmp < start+1;
% idx2 = tmp >= start+1 & tmp < start+2;
% % idx2 = idx;
% tmp2 = find(idx);
% offset = -vecPOS_F(tmp2(1));
% hold on
% plot(vecPOS_F(idx) + offset, (CT_fixed(idx) + CT_fixed(idx2))./2, '-.m');
% hold off

legend('Experimental Data','DDE Method','Location','NorthWest','FontSize',8)
xlabel('Azimuth Location (Degrees)');
ylabel('Thrust Coefficient');
% title('\mu = 0.3, \alpha = 0, RPM = 3000, \DeltaT = 0.0005')
WH = [4.5*2 5];
fcnFIG2LATEX(hFig2, 'tmotor_time_2.pdf', WH)

%% FFT

% OPERA
start = 41;
t_op = [start:valMAXTIME].*valDELTIME;
x_op = CT_U(start:end);
n = length(x_op);
y_op = fft(detrend(x_op));
f_op = (1/valDELTIME).*(0:length(y_op)-1)/length(y_op);
p_op = abs(y_op)/n;

% Tunnel
CT_tunnel = lbf_N.*FT(:,3);
CT_tunnel = CT_tunnel./(rho.*(pi.*((valDIAM/2).^2)).*(((valDIAM/2).*(valRPM.*(pi/30))).^2));

rate = 10000; % 10k Hz
n = length(CT_tunnel);
y_t = fft(detrend(CT_tunnel), n);
f_t = rate.*(0:length(y_t)-1)/(length(y_t));
p_t = abs(y_t)/n;

hFig68 = figure(68);
clf(68);
hold on
plot(f_t, p_t)
% yyaxis right
scatter(f_op, p_op);
hold off

xlim([0 1000])
xlabel('Hz')
ylabel('Signal Power');

grid minor
box on
% axis tight

legend('RU Test Data', 'DDE Method')

%%
% hFig20 = figure(20);
% clf(20);
% 
% load('Alpha 0 Results/TMotor_Relaxed_J0.3_0.0005.mat', 'matDVELIFT_DIR', 'matDVEDRAG_DIR','matDVESIDE_DIR','matLIFTFREE','matLIFTIND','matSIDEFREE','matSIDEIND','matDRAGIND','valDIAM','valRPM','valDENSITY', 'valUINF', 'valAREA', 'matINTCIRC', 'strATYPE');
% 
% cd ../../
% [CT_U, CL_U, matDGAMMADT] = fcnDGAMMADT(1, valDELTIME, strATYPE, matINTCIRC, valDENSITY, valRPM, valDIAM, valAREA, valUINF, matLIFTFREE, matLIFTIND, matDRAGIND, matSIDEFREE, matSIDEIND, matDVELIFT_DIR, matDVEDRAG_DIR, matDVESIDE_DIR);
% cd("Stuff/TMotor Study/")
% 
% matDVETHRUST = dot(matDVELIFT_DIR.*permute(matLIFTFREE + matLIFTIND + matDGAMMADT, [2 3 1]), repmat([0 0 1], 20, 1, valMAXTIME), 2) ...
%     + dot(matDVEDRAG_DIR.*permute(matDRAGIND, [2 3 1]), repmat([0 0 1], 20, 1, valMAXTIME), 2) ...
%     + dot(matDVESIDE_DIR.*permute(matSIDEFREE + matSIDEIND, [2 3 1]), repmat([0 0 1], 20, 1, valMAXTIME), 2);
% matDVETHRUST = matDVETHRUST./(valDENSITY.*(pi.*((valDIAM/2).^2)).*(((valDIAM/2).*(valRPM.*(pi/30))).^2));
% 
% WH = [4.5 5];
% idx3 = vecDVESURFACE(vecTEDVE);
% tmp5 = vecPOS_R(idx) + offset;
% idx3 = (vecDVESURFACE(vecTEDVE));
% 
% hold on
% plot(tmp5, reshape((sum(matDVETHRUST(idx3==2,:,idx),1) + sum(matDVETHRUST(idx3==2,:,idx2),1))./2,[],1,1), '--b')
% plot(tmp5, reshape((sum(matDVETHRUST(idx3==1,:,idx),1) + sum(matDVETHRUST(idx3==1,:,idx2),1))./2,[],1,1), '-.r')
% plot(tmp5, reshape((sum(matDVETHRUST(:,:,idx),1) + sum(matDVETHRUST(:,:,idx2),1))./2,[],1,1), '-k')
% hold off
% axis tight
% grid minor
% box on
% 
% xlabel('Position (Degrees)');
% ylabel('Thrust (N)');
% legend('Blade A','Blade B', 'Combined','Location','NorthWest','FontSize',10)
% 
% WH = [4.5*2 5];
% fcnFIG2LATEX(hFig20, 'blade_thrust_relaxed1.pdf', WH)

hFig20 = figure(20);
clf(20);

load('Alpha 0 Results/TMotor_Relaxed_J0.3_0.0005.mat', 'matDVELIFT_DIR', 'matDVEDRAG_DIR','matDVESIDE_DIR','matLIFTFREE','matLIFTIND','matSIDEFREE','matSIDEIND','matDRAGIND','valDIAM','valRPM','valDENSITY', 'valUINF', 'valAREA', 'matINTCIRC', 'strATYPE');

cd ../../
[CT_U, CL_U, matDGAMMADT] = fcnDGAMMADT(1, valDELTIME, strATYPE, matINTCIRC, valDENSITY, valRPM, valDIAM, valAREA, valUINF, matLIFTFREE, matLIFTIND, matDRAGIND, matSIDEFREE, matSIDEIND, matDVELIFT_DIR, matDVEDRAG_DIR, matDVESIDE_DIR);
cd("Stuff/TMotor Study/")

subplot(3,1,2)
matDVETHRUST = dot(matDVELIFT_DIR.*permute(matLIFTIND, [2 3 1]), repmat([0 0 1], 20, 1, valMAXTIME), 2) ...
    + dot(matDVEDRAG_DIR.*permute(matDRAGIND, [2 3 1]), repmat([0 0 1], 20, 1, valMAXTIME), 2) ...
    + dot(matDVESIDE_DIR.*permute(matSIDEIND, [2 3 1]), repmat([0 0 1], 20, 1, valMAXTIME), 2);
matDVETHRUST = matDVETHRUST./(valDENSITY.*(pi.*((valDIAM/2).^2)).*(((valDIAM/2).*(valRPM.*(pi/30))).^2));

WH = [4.5 5];
idx3 = vecDVESURFACE(vecTEDVE);
tmp5 = vecPOS_R(idx) + offset;
idx3 = (vecDVESURFACE(vecTEDVE));

hold on
plot(tmp5, reshape((sum(matDVETHRUST(idx3==2,:,idx),1) + sum(matDVETHRUST(idx3==2,:,idx2),1))./2,[],1,1), '--b')
plot(tmp5, reshape((sum(matDVETHRUST(idx3==1,:,idx),1) + sum(matDVETHRUST(idx3==1,:,idx2),1))./2,[],1,1), '-.r')
plot(tmp5, reshape((sum(matDVETHRUST(:,:,idx),1) + sum(matDVETHRUST(:,:,idx2),1))./2,[],1,1), '-k')
hold off
axis tight
grid minor
box on

title('Induced');
xlabel('Position (Degrees)');
ylabel('Thrust (N)');
legend('Blade A','Blade B', 'Combined','Location','NorthWest','FontSize',10)

subplot(3,1,1)
matDVETHRUST = dot(matDVELIFT_DIR.*permute(matLIFTFREE, [2 3 1]), repmat([0 0 1], 20, 1, valMAXTIME), 2) ...
    + dot(matDVESIDE_DIR.*permute(matSIDEFREE, [2 3 1]), repmat([0 0 1], 20, 1, valMAXTIME), 2);
matDVETHRUST = matDVETHRUST./(valDENSITY.*(pi.*((valDIAM/2).^2)).*(((valDIAM/2).*(valRPM.*(pi/30))).^2));

WH = [4.5 5];
idx3 = vecDVESURFACE(vecTEDVE);
tmp5 = vecPOS_R(idx) + offset;
idx3 = (vecDVESURFACE(vecTEDVE));

hold on
plot(tmp5, reshape((sum(matDVETHRUST(idx3==2,:,idx),1) + sum(matDVETHRUST(idx3==2,:,idx2),1))./2,[],1,1), '--b')
plot(tmp5, reshape((sum(matDVETHRUST(idx3==1,:,idx),1) + sum(matDVETHRUST(idx3==1,:,idx2),1))./2,[],1,1), '-.r')
plot(tmp5, reshape((sum(matDVETHRUST(:,:,idx),1) + sum(matDVETHRUST(:,:,idx2),1))./2,[],1,1), '-k')
hold off
axis tight
grid minor
box on

title('Freestream')
xlabel('Position (Degrees)');
ylabel('Thrust (N)');
legend('Blade A','Blade B', 'Combined','Location','NorthWest','FontSize',10)

subplot(3,1,3)
matDVETHRUST = dot(matDVELIFT_DIR.*permute(matDGAMMADT, [2 3 1]), repmat([0 0 1], 20, 1, valMAXTIME), 2);
matDVETHRUST = matDVETHRUST./(valDENSITY.*(pi.*((valDIAM/2).^2)).*(((valDIAM/2).*(valRPM.*(pi/30))).^2));

WH = [4.5 5];
idx3 = vecDVESURFACE(vecTEDVE);
tmp5 = vecPOS_R(idx) + offset;
idx3 = (vecDVESURFACE(vecTEDVE));

title('Unsteady');
hold on
plot(tmp5, reshape((sum(matDVETHRUST(idx3==2,:,idx),1) + sum(matDVETHRUST(idx3==2,:,idx2),1))./2,[],1,1), '--b')
plot(tmp5, reshape((sum(matDVETHRUST(idx3==1,:,idx),1) + sum(matDVETHRUST(idx3==1,:,idx2),1))./2,[],1,1), '-.r')
plot(tmp5, reshape((sum(matDVETHRUST(:,:,idx),1) + sum(matDVETHRUST(:,:,idx2),1))./2,[],1,1), '-k')
hold off
axis tight
grid minor
box on

xlabel('Position (Degrees)');
ylabel('Thrust (N)');
legend('Blade A','Blade B', 'Combined','Location','NorthWest','FontSize',10)




