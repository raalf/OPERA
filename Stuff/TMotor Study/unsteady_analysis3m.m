clc
clear

load('Tunnel Testing\2020-01-30\30-Jan-2020 16.15.50_Scorpion_KDE_T-Motor 18in_RPM3000_Alpha0_20.7492.mat', ... % Alpha 0, 0.2889
    'Angle', 'lbf_N', 'FT', 'rho', 'valDIAM', 'valRPM', 'vecPOS_TUNNEL_OG', 'dataRate');

% load('Tunnel Testing\2020-01-30\30-Jan-2020 16.40.19_Scorpion_KDE_T-Motor 18in_RPM3000_Alpha15_7.9674.mat', ... % Alpha 15, 0.1109
%     'Angle', 'lbf_N', 'FT', 'rho', 'valDIAM', 'valRPM', 'vecPOS_TUNNEL_OG', 'dataRate');

CT_tunnel_raw = lbf_N.*FT(:,3);
CT_tunnel_raw = CT_tunnel_raw./(rho.*(pi.*((valDIAM/2).^2)).*(((valDIAM/2).*(valRPM.*(pi/30))).^2));


%% Getting OPERA data
load('Alpha 0 Results/TMotor_Relaxed_J0.3_0.0005_440.mat', 'CT_U', 'CT', 'valDELTIME', 'matDVE', 'matVLST', 'vecHUB', 'vecDVESURFACE', 'matDGAMMADT', 'matINTCIRC', 'vecTEDVE', 'valMAXTIME', 'matSPANDIR', 'valRPM', 'vecDVETHRUST')
% load('Alpha 0 Results/TMotor_Fixed_J0.3_0.0005_440.mat', 'CT_U', 'CT', 'valDELTIME', 'matDVE', 'matVLST', 'vecHUB', 'vecDVESURFACE', 'matDGAMMADT', 'matINTCIRC', 'vecTEDVE', 'valMAXTIME', 'matSPANDIR', 'valRPM', 'vecDVETHRUST')

% load('Alpha 15 Results/TMotor_Relaxed_J0.1109_0.0005_440.mat', 'CT_U', 'CT', 'valDELTIME', 'matDVE', 'matVLST', 'vecHUB', 'vecDVESURFACE', 'matDGAMMADT', 'matINTCIRC', 'vecTEDVE', 'valMAXTIME', 'matSPANDIR', 'valRPM', 'vecDVETHRUST')
% load('Alpha 15 Results/TMotor_Fixed_J0.1109_0.0005_440.mat', 'CT_U', 'CT', 'valDELTIME', 'matDVE', 'matVLST', 'vecHUB', 'vecDVESURFACE', 'matDGAMMADT', 'matINTCIRC', 'vecTEDVE', 'valMAXTIME', 'matSPANDIR', 'valRPM', 'vecDVETHRUST')

CT_relaxed = CT_U(~isnan(CT_U));

%% Plotting
hFig2 = figure(2);
clf(2);

binAng = linspace(0, 360, 30);
for i = 1:length(binAng)
    rng = 0.25;
    idx = vecPOS_TUNNEL_OG >= binAng(i) - rng & vecPOS_TUNNEL_OG <= binAng(i) + rng;
    
    binAvg(i) = mean(CT_tunnel_raw(idx));
    binMax(i) = max(CT_tunnel_raw(idx));
    binMin(i) = min(CT_tunnel_raw(idx));
end
plot(binAng, binAvg, ':sk');

%---------------------------------------------

Fp = 1/dataRate;
Fst = 1000/dataRate;
d = fdesign.lowpass('N,Fp,Fst',3,Fp,Fst);
Hd = design(d);
CT_tunnel2 = filtfilt(Hd.Numerator,1,detrend(CT_tunnel_raw)) + mean(CT_tunnel_raw);

binAng = linspace(0, 360, 30);
binAvg = [];
binMax = [];
binMin = [];
for i = 1:length(binAng)
    rng = 0.25;
    idx = vecPOS_TUNNEL_OG >= binAng(i) - rng & vecPOS_TUNNEL_OG <= binAng(i) + rng;
    binAvg(i) = mean(CT_tunnel2(idx));
    binMax(i) = max(CT_tunnel2(idx));
    binMin(i) = min(CT_tunnel2(idx));
end
hold on
plot(binAng, binAvg, '-.^m');
hold off

title('SHIFTED Alpha 15, \mu = 0.11')

%%
deg_per_ts = valRPM.*(pi/30).*(180/pi).*valDELTIME;
tmp_offset = 0;
vecPOS_R = [tmp_offset:(length(CT_relaxed)-1 + tmp_offset)]'.*deg_per_ts + 90;
vecPOS = mod(vecPOS_R,360);

binAng = linspace(0, 360, 36);
binAvg = [];
binMax = [];
binMin = [];
CT_U(1:80) = nan;
for i = 1:length(binAng)
    rng = 4.5;
    idx = vecPOS >= binAng(i) - rng & vecPOS <= binAng(i) + rng;
    if any(idx)
        binAvg(i) = nanmean(CT_U(idx));
        binMax(i) = nanmax(CT_U(idx));
        binMin(i) = nanmin(CT_U(idx));
    else
        binAvg(i) = nan;
        binMax(i) = nan;
        binMin(i) = nan;       
    end
end
hold on
plot(binAng, binAvg, '--b');
hold off

% tmp = vecPOS_R./360;
% start = 9;
% 
% idx = tmp >= start & tmp < start+1;
% idx2 = tmp >= start+1 & tmp < start+2;
% % idx2 = idx;
% CT_relaxed2 = CT(~isnan(CT));
% tmp2 = find(idx);
% offset = -vecPOS_R(tmp2(1));
% hold on
% plot(vecPOS_R(idx) + offset, (CT_relaxed(idx) + CT_relaxed(idx2))./2, '--b');
% % plot(vecPOS_R(idx) + offset, (CT_relaxed2(idx) + CT_relaxed2(idx2))./2, '-.r');

% hold off
grid minor
box on
axis tight

legend('Experimental Data','Experimental Data (0 - 2000 Hz)','DDE Method','Location','NorthEast')
xlabel('Azimuth Location (Degrees)');
ylabel('Thrust Coefficient');
% title('\mu = 0.3, \alpha = 0, RPM = 3000, \DeltaT = 0.0005')

WH = [4.5 5];
% fcnFIG2LATEX(hFig2, 'tmotor_time_15.pdf', WH)


