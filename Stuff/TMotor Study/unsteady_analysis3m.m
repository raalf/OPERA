clc
clear

% load('Tunnel Testing\2019-10-31\31-Oct-2019 12.41.35_Scorpion_ASI_T-Motor 18in_RPM3000_Alpha30_21.5736.mat') % Alpha 0, 0.3004
load('Tunnel Testing\2019-11-06\06-Nov-2019 13.02.45_Scorpion_ASI_T-Motor 18in_RPM3000_Alpha15_8.0104.mat') % Alpha 15, 0.1115

CT_tunnel_raw = lbf_N.*FT(:,3);
CT_tunnel_raw = CT_tunnel_raw./(rho.*(pi.*((valDIAM/2).^2)).*(((valDIAM/2).*(valRPM.*(pi/30))).^2));

% close all

%% Getting OPERA data
% load('Alpha 0 Results/TMotor_Relaxed_J0.3_0.0005.mat', 'CT_U', 'CT', 'valDELTIME', 'matDVE', 'matVLST', 'vecHUB', 'vecDVESURFACE', 'matDGAMMADT', 'matINTCIRC', 'vecTEDVE', 'valMAXTIME', 'matSPANDIR', 'valRPM', 'vecDVETHRUST')
load('Alpha 15 Results/TMotor_Relaxed_J0.1115_0.0005.mat', 'CT_U', 'CT', 'valDELTIME', 'matDVE', 'matVLST', 'vecHUB', 'vecDVESURFACE', 'matDGAMMADT', 'matINTCIRC', 'vecTEDVE', 'valMAXTIME', 'matSPANDIR', 'valRPM', 'vecDVETHRUST')

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
% errorbar(binAng, binAvg, binMin - binAvg, binMax - binAvg, ':sk');
plot(binAng, binAvg, ':sk');

%---------------------------------------------

Fp = 1/dataRate;
Fst = 3000/dataRate;
d = fdesign.lowpass('N,Fp,Fst',3,Fp,Fst);
Hd = design(d);
CT_tunnel2 = filtfilt(Hd.Numerator,1,detrend(CT_tunnel_raw)) + mean(CT_tunnel_raw);

% CT_tunnel2 = lowpass(CT_tunnel_raw, 300, dataRate, 'Steepness', 0.999, 'StopbandAttenuation', 10);
% CT_tunnel2 = lowpass(CT_tunnel_raw, 300, dataRate)
[~,idx3] = sort(Angle,'ascend');
CT_tunnel2 = CT_tunnel2(idx3);

binAng = linspace(0, 360, 30);
for i = 1:length(binAng)
    rng = 0.25;
    idx = vecPOS_TUNNEL_OG >= binAng(i) - rng & vecPOS_TUNNEL_OG <= binAng(i) + rng;
    binAvg(i) = mean(CT_tunnel2(idx));
    binMax(i) = max(CT_tunnel2(idx));
    binMin(i) = min(CT_tunnel2(idx));
end
hold on
% errorbar(binAng, binAvg, binMin - binAvg, binMax - binAvg, '-.^m');
plot(binAng, binAvg, '-.^m');
hold off

%%

tmp = ((vecPOS_R + 90)./360);
% start = 2;
start = 1.9;
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

% legend('Experimental Data','DDE Method','Location','South','FontSize',8)
legend('Experimental Data','Experimental Data (0 - 500 Hz)','DDE Method','Location','SouthWest')
xlabel('Azimuth Location (Degrees)');
ylabel('Thrust Coefficient');
% title('\mu = 0.3, \alpha = 0, RPM = 3000, \DeltaT = 0.0005')
WH = [4.5 5];
% fcnFIG2LATEX(hFig2, 'tmotor_time_15.pdf', WH)

% title('Alpha 15, J = 0.115')

%%
% hFig20 = figure(20);
% clf(20);
% 
% load('Alpha 0 Results/TMotor_Relaxed_J0.3_0.0005.mat', 'matDVELIFT_DIR', 'matDVEDRAG_DIR','matDVESIDE_DIR','matLIFTFREE','matLIFTIND','matSIDEFREE','matSIDEIND','matDRAGIND','valDIAM','valRPM','valDENSITY', 'valUINF', 'valAREA', 'matINTCIRC', 'strATYPE', 'matSPANDIR_ALL', 'matUINF_ALL', 'vecTE', 'vecTEDVE');
% % load('Alpha 15 Results/TMotor_Relaxed_J0.1115_0.0005.mat', 'matDVELIFT_DIR', 'matDVEDRAG_DIR','matDVESIDE_DIR','matLIFTFREE','matLIFTIND','matSIDEFREE','matSIDEIND','matDRAGIND','valDIAM','valRPM','valDENSITY', 'valUINF', 'valAREA', 'matINTCIRC', 'strATYPE', 'matSPANDIR_ALL', 'matUINF_ALL', 'vecTE', 'vecTEDVE');
% 
% cd ../../
% [CT_U, CL_U, matDGAMMADT] = fcnDGAMMADT(1, valDELTIME, strATYPE, matINTCIRC, valDENSITY, valRPM, valDIAM, valAREA, valUINF, matLIFTFREE, matLIFTIND, matDRAGIND, matSIDEFREE, matSIDEIND, matDVELIFT_DIR, matDVEDRAG_DIR, matDVESIDE_DIR, matSPANDIR_ALL, matUINF_ALL, vecTE, vecTEDVE);
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
% plot(tmp5, reshape((sum(matDVETHRUST(:,:,idx),1) + sum(matDVETHRUST(:,:,idx2),1))./2,[],1,1), '-k')
% plot(tmp5, reshape((sum(matDVETHRUST(idx3==2,:,idx),1) + sum(matDVETHRUST(idx3==2,:,idx2),1))./2,[],1,1), '--b')
% plot(tmp5, reshape((sum(matDVETHRUST(idx3==1,:,idx),1) + sum(matDVETHRUST(idx3==1,:,idx2),1))./2,[],1,1), '-.r')
% hold off
% axis tight
% grid minor
% box on
% 
% xlabel('Azimuth Location (Degrees)');
% ylabel('Thrust Coefficient');
% % legend('Combined', 'Blade A','Blade B','Location','West','FontSize',8)
% 
% WH = [4.5 5];
% % fcnFIG2LATEX(hFig20, 'blade_thrust_relaxed_15.pdf', WH)
% fcnFIG2LATEX(hFig20, 'blade_thrust_relaxed_0.pdf', WH)
% 
% 
% % hFig20 = figure(20);
% % clf(20);
% % 
% % % load('Alpha 0 Results/TMotor_Relaxed_J0.3_0.0005.mat', 'matDVELIFT_DIR', 'matDVEDRAG_DIR','matDVESIDE_DIR','matLIFTFREE','matLIFTIND','matSIDEFREE','matSIDEIND','matDRAGIND','valDIAM','valRPM','valDENSITY', 'valUINF', 'valAREA', 'matINTCIRC', 'strATYPE', 'matSPANDIR_ALL', 'matUINF_ALL', 'vecTE', 'vecTEDVE');
% % % load('Alpha 0 Results/TMotor_Relaxed_J0.1346_0.0005.mat', 'matDVELIFT_DIR', 'matDVEDRAG_DIR','matDVESIDE_DIR','matLIFTFREE','matLIFTIND','matSIDEFREE','matSIDEIND','matDRAGIND','valDIAM','valRPM','valDENSITY', 'valUINF', 'valAREA', 'matINTCIRC', 'strATYPE', 'matSPANDIR_ALL', 'matUINF_ALL', 'vecTE', 'vecTEDVE');
% % 
% % cd ../../
% % [CT_U, CL_U, matDGAMMADT] = fcnDGAMMADT(1, valDELTIME, strATYPE, matINTCIRC, valDENSITY, valRPM, valDIAM, valAREA, valUINF, matLIFTFREE, matLIFTIND, matDRAGIND, matSIDEFREE, matSIDEIND, matDVELIFT_DIR, matDVEDRAG_DIR, matDVESIDE_DIR, matSPANDIR_ALL, matUINF_ALL, vecTE, vecTEDVE);
% % cd("Stuff/TMotor Study/")
% % 
% % subplot(3,1,2)
% % matDVETHRUST = dot(matDVELIFT_DIR.*permute(matLIFTIND, [2 3 1]), repmat([0 0 1], 20, 1, valMAXTIME), 2) ...
% %     + dot(matDVEDRAG_DIR.*permute(matDRAGIND, [2 3 1]), repmat([0 0 1], 20, 1, valMAXTIME), 2) ...
% %     + dot(matDVESIDE_DIR.*permute(matSIDEIND, [2 3 1]), repmat([0 0 1], 20, 1, valMAXTIME), 2);
% % matDVETHRUST = matDVETHRUST./(valDENSITY.*(pi.*((valDIAM/2).^2)).*(((valDIAM/2).*(valRPM.*(pi/30))).^2));
% % 
% % WH = [4.5 5];
% % idx3 = vecDVESURFACE(vecTEDVE);
% % tmp5 = vecPOS_R(idx) + offset;
% % idx3 = (vecDVESURFACE(vecTEDVE));
% % 
% % hold on
% % plot(tmp5, reshape((sum(matDVETHRUST(idx3==2,:,idx),1) + sum(matDVETHRUST(idx3==2,:,idx2),1))./2,[],1,1), '--b')
% % plot(tmp5, reshape((sum(matDVETHRUST(idx3==1,:,idx),1) + sum(matDVETHRUST(idx3==1,:,idx2),1))./2,[],1,1), '-.r')
% % plot(tmp5, reshape((sum(matDVETHRUST(:,:,idx),1) + sum(matDVETHRUST(:,:,idx2),1))./2,[],1,1), '-k')
% % hold off
% % axis tight
% % grid minor
% % box on
% % 
% % title('Induced');
% % xlabel('Position (Degrees)');
% % ylabel('Thrust (N)');
% % legend('Blade A','Blade B', 'Combined','Location','NorthWest','FontSize',10)
% % 
% % subplot(3,1,1)
% % matDVETHRUST = dot(matDVELIFT_DIR.*permute(matLIFTFREE, [2 3 1]), repmat([0 0 1], 20, 1, valMAXTIME), 2) ...
% %     + dot(matDVESIDE_DIR.*permute(matSIDEFREE, [2 3 1]), repmat([0 0 1], 20, 1, valMAXTIME), 2);
% % matDVETHRUST = matDVETHRUST./(valDENSITY.*(pi.*((valDIAM/2).^2)).*(((valDIAM/2).*(valRPM.*(pi/30))).^2));
% % 
% % WH = [4.5 5];
% % idx3 = vecDVESURFACE(vecTEDVE);
% % tmp5 = vecPOS_R(idx) + offset;
% % idx3 = (vecDVESURFACE(vecTEDVE));
% % 
% % hold on
% % plot(tmp5, reshape((sum(matDVETHRUST(idx3==2,:,idx),1) + sum(matDVETHRUST(idx3==2,:,idx2),1))./2,[],1,1), '--b')
% % plot(tmp5, reshape((sum(matDVETHRUST(idx3==1,:,idx),1) + sum(matDVETHRUST(idx3==1,:,idx2),1))./2,[],1,1), '-.r')
% % plot(tmp5, reshape((sum(matDVETHRUST(:,:,idx),1) + sum(matDVETHRUST(:,:,idx2),1))./2,[],1,1), '-k')
% % hold off
% % axis tight
% % grid minor
% % box on
% % 
% % title('Freestream')
% % xlabel('Position (Degrees)');
% % ylabel('Thrust (N)');
% % legend('Blade A','Blade B', 'Combined','Location','NorthWest','FontSize',10)
% % 
% % subplot(3,1,3)
% % matDVETHRUST = dot(matDVELIFT_DIR.*permute(matDGAMMADT, [2 3 1]), repmat([0 0 1], 20, 1, valMAXTIME), 2);
% % matDVETHRUST = matDVETHRUST./(valDENSITY.*(pi.*((valDIAM/2).^2)).*(((valDIAM/2).*(valRPM.*(pi/30))).^2));
% % 
% % WH = [4.5 5];
% % idx3 = vecDVESURFACE(vecTEDVE);
% % tmp5 = vecPOS_R(idx) + offset;
% % idx3 = (vecDVESURFACE(vecTEDVE));
% % 
% % title('Unsteady');
% % hold on
% % plot(tmp5, reshape((sum(matDVETHRUST(idx3==2,:,idx),1) + sum(matDVETHRUST(idx3==2,:,idx2),1))./2,[],1,1), '--b')
% % plot(tmp5, reshape((sum(matDVETHRUST(idx3==1,:,idx),1) + sum(matDVETHRUST(idx3==1,:,idx2),1))./2,[],1,1), '-.r')
% % plot(tmp5, reshape((sum(matDVETHRUST(:,:,idx),1) + sum(matDVETHRUST(:,:,idx2),1))./2,[],1,1), '-k')
% % hold off
% % axis tight
% % grid minor
% % box on
% % 
% % xlabel('Position (Degrees)');
% % ylabel('Thrust (N)');
% % legend('Blade A','Blade B', 'Combined','Location','NorthWest','FontSize',10)
% % 
% % 
% % %%



