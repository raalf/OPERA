clc
clear

load('Alpha 0 Results/TMotor_Relaxed_J0.3_0.0005_440.mat', 'CT_U', 'CT', 'valDELTIME', 'matDVE', 'matVLST', 'vecHUB', 'vecDVESURFACE', 'matDGAMMADT', 'matINTCIRC', 'vecTEDVE', 'valMAXTIME', 'matSPANDIR', 'valRPM', 'vecDVETHRUST');
load('G:\GIT\opera\Stuff\TMotor Study\Tunnel Testing\2019-10-31\31-Oct-2019 12.41.35_Scorpion_ASI_T-Motor 18in_RPM3000_Alpha30_21.5736.mat') % 0.3004
CT_relaxed = CT_U(~isnan(CT_U));
% CT_relaxed = CT(~isnan(CT));

deg_per_ts = valRPM.*(pi/30).*(180/pi).*valDELTIME;
tmp_offset = 0;
vecPOS_R = [tmp_offset:(length(CT_relaxed)-1 + tmp_offset)]'.*deg_per_ts;

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