clc
clear

% load('Alpha 0 Results/TMotor_Relaxed_J0.3_0.0005_440.mat', 'CT_U', 'CT', 'valDELTIME', 'matDVE', 'matVLST', 'vecHUB', 'vecDVESURFACE', 'matDGAMMADT', 'matINTCIRC', 'vecTEDVE', 'valMAXTIME', 'matSPANDIR', 'valRPM', 'vecDVETHRUST');
% load('Tunnel Testing\2019-10-31\31-Oct-2019 12.41.35_Scorpion_ASI_T-Motor 18in_RPM3000_Alpha30_21.5736.mat') % 0.3004

load('Alpha 15 Results/TMotor_Relaxed_J0.1115_0.0005_440.mat', 'CT_U', 'CT', 'valDELTIME', 'matDVE', 'matVLST', 'vecHUB', 'vecDVESURFACE', 'matDGAMMADT', 'matINTCIRC', 'vecTEDVE', 'valMAXTIME', 'matSPANDIR', 'valRPM', 'vecDVETHRUST');
load('Tunnel Testing\2019-11-06\06-Nov-2019 13.02.45_Scorpion_ASI_T-Motor 18in_RPM3000_Alpha15_8.0104.mat') % Alpha 15, 0.1115

close all
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

%% Julia Cole's method
% bpFilt = designfilt('bandpassfir','FilterOrder',20,'CutoffFrequency1',10,'CutoffFrequency2',200,'SampleRate',dataRate);
bpFilt = designfilt('bandpassiir','FilterOrder',20,'HalfPowerFrequency1',75,'HalfPowerFrequency2',250,'SampleRate',dataRate);
% bpFilt = designfilt('lowpassiir','FilterOrder',20,'PassbandFrequency',400,'Passbandripple',0.1,'SampleRate',dataRate);
CT_tunnel2 = filtfilt(bpFilt,CT_tunnel);


% [B,A] = butter(order, frequency, 'low');
% Filterddata = filter(B,A,datatobefiltered);

% cutoff = 400/(dataRate/2);
% [B,A] = butter(3, cutoff, 'low');
% CT_tunnel2 = filter(B,A,CT_tunnel);

%%


% rate = 10000; % 10k Hz
% n = length(CT_tunnel);
% y_t = fft(detrend(CT_tunnel), n);
% f_t = rate.*(0:length(y_t)-1)/(length(y_t));
% p_t = abs(y_t)/n;

hFig68 = figure(68);
clf(68);
hold on
% plot(f_t, p_t, '-b')
% yyaxis right
plot(f_op, p_op, 'sk');

hold off


%------
rate = 10000; % 10k Hz
n = length(CT_tunnel);
y_t = fft(detrend(CT_tunnel2), n);
f_t = rate.*(0:length(y_t)-1)/(length(y_t));
p_t = abs(y_t)/n;

hold on
plot(f_t, p_t, '--k')
hold off
xlim([0 400])

%--------

% xlim([0 1000])
xlabel('Hz')
ylabel('Signal Power');

grid minor
box on
% axis tight

% legend('Experimental', 'DDE Method')

WH = [4.5 5];
% fcnFIG2LATEX(hFig68, 'tmotor_fft_0.pdf', WH)

%%
% clc
% clear
% 
% load('Alpha 15 Results/TMotor_Relaxed_J0.1115_0.0005_440.mat', 'CT_U', 'CT', 'valDELTIME', 'matDVE', 'matVLST', 'vecHUB', 'vecDVESURFACE', 'matDGAMMADT', 'matINTCIRC', 'vecTEDVE', 'valMAXTIME', 'matSPANDIR', 'valRPM', 'vecDVETHRUST');
% load('G:\GIT\opera\Stuff\TMotor Study\Tunnel Testing\2019-11-06\06-Nov-2019 13.02.45_Scorpion_ASI_T-Motor 18in_RPM3000_Alpha15_8.0104.mat') % Alpha 15, 0.1115
% 
% close all
% CT_relaxed = CT_U(~isnan(CT_U));
% % CT_relaxed = CT(~isnan(CT));
% 
% deg_per_ts = valRPM.*(pi/30).*(180/pi).*valDELTIME;
% tmp_offset = 0;
% vecPOS_R = [tmp_offset:(length(CT_relaxed)-1 + tmp_offset)]'.*deg_per_ts;
% 
% % OPERA
% start = 41;
% t_op = [start:valMAXTIME].*valDELTIME;
% x_op = CT_U(start:end);
% n = length(x_op);
% y_op = fft(detrend(x_op));
% f_op = (1/valDELTIME).*(0:length(y_op)-1)/length(y_op);
% p_op = abs(y_op)/n;
% 
% % Tunnel
% CT_tunnel = lbf_N.*FT(:,3);
% CT_tunnel = CT_tunnel./(rho.*(pi.*((valDIAM/2).^2)).*(((valDIAM/2).*(valRPM.*(pi/30))).^2));
% 
% rate = 10000; % 10k Hz
% n = length(CT_tunnel);
% y_t = fft(detrend(CT_tunnel), n);
% f_t = rate.*(0:length(y_t)-1)/(length(y_t));
% p_t = abs(y_t)/n;
% 
% hFig67 = figure(67);
% clf(67);
% hold on
% plot(f_t, p_t)
% % yyaxis right
% plot(f_op, p_op, 'sk');
% hold off
% 
% xlim([0 1000])
% xlabel('Hz')
% ylabel('Signal Power');
% 
% grid minor
% box on
% % axis tight
% 
% legend('Experimental', 'DDE Method')
% 
% WH = [4.5 5];
% fcnFIG2LATEX(hFig67, 'tmotor_fft_15.pdf', WH)
