clc
clear

%% Tunnel
load('Tunnel Testing\2020-01-30\30-Jan-2020 16.15.50_Scorpion_KDE_T-Motor 18in_RPM3000_Alpha0_20.7492.mat', ... % Alpha 0, 0.2889
'Angle', 'lbf_N', 'FT', 'rho', 'valDIAM', 'valRPM', 'vecPOS_TUNNEL_OG', 'dataRate');

fs = dataRate; % sampling rate
rps = valRPM/60; % rotation rate
T = 4*60; % total time
N = fs*T; % total samples
t = linspace(0, T, N);

FFTParam.FFTLength = fs; % notice how your BR peaks change when you adjust fft length here, e.g. fs*10
FFTParam.WinType = 'hanning'; % play with flattopwin or rectwin to see difference windowing makes on spectra
FFTParam.SampRate = fs; 
FFTParam.DFlag = 'none'; 
FFTParam.PercOverlap = 50;

% Alpha 0
CT_tunnel = lbf_N.*FT(:,3);
CT_tunnel = CT_tunnel./(rho.*(pi.*((valDIAM/2).^2)).*(((valDIAM/2).*(valRPM.*(pi/30))).^2));
InputData = detrend(CT_tunnel);
[Pxy,COH, F, Navg]=CalculatePSDandCohs(InputData, FFTParam);

hFig1 = figure(1);
clf(1);
plot(F, 10*log10(Pxy), '--k')

% Alpha 15
load('Tunnel Testing\2020-01-30\30-Jan-2020 16.40.19_Scorpion_KDE_T-Motor 18in_RPM3000_Alpha15_7.9674.mat', ... % Alpha 15, 0.1109
    'Angle', 'lbf_N', 'FT', 'rho', 'valDIAM', 'valRPM', 'vecPOS_TUNNEL_OG', 'dataRate');
CT_tunnel = lbf_N.*FT(:,3);
CT_tunnel = CT_tunnel./(rho.*(pi.*((valDIAM/2).^2)).*(((valDIAM/2).*(valRPM.*(pi/30))).^2));
InputData = detrend(CT_tunnel);
[Pxy,COH, F, Navg]=CalculatePSDandCohs(InputData, FFTParam);

hFig2 = figure(2);
clf(2);
plot(F, 10*log10(Pxy), '--k')

%% DDE
% Alpha 0
load('G:\GIT\opera\Stuff\TMotor Study\Alpha 0 Results\TMotor_Relaxed_J0.3_0.0005_440.mat')

deg_per_ts = valRPM.*(pi/30).*(180/pi).*valDELTIME;
vecPOS_R = [0:(length(CT_U) - 1)]'.*deg_per_ts;
vecPOS = mod(vecPOS_R,360);

npts = 40;
cutoff = 2;
n = ceil((4*60)*((valRPM/60)./(11 - cutoff))); % 4 minutes?
CT_ext = repmat(CT_U((cutoff*npts + 1):end),n,1);

fs = 1/valDELTIME; % sampling rate
rps = valRPM/60; % rotation rate
T = length(CT_ext).*valDELTIME; % total time
N = fs*T; % total samples
t = linspace(0, T, N);

FFTParam.FFTLength = fs; % notice how your BR peaks change when you adjust fft length here, e.g. fs*10
FFTParam.WinType = 'rectwin'; % play with flattopwin or rectwin to see difference windowing makes on spectra
FFTParam.SampRate = fs; 
FFTParam.DFlag = 'none'; 
FFTParam.PercOverlap = 50;

InputData = detrend(CT_ext);
[Pxy,COH, F, Navg] = CalculatePSDandCohs(InputData, FFTParam);

figure(1);
hold on
plot(F, 10*log10(Pxy), '-.b')
hold off

xlabel('Frequency (Hz)')
ylabel('Power/Frequency (dB/Hz)')
grid minor
box on
axis tight
title('Alpha 0, \mu = 0.3')

xlim([0 1000])

% Alpha 15
load('G:\GIT\opera\Stuff\TMotor Study\Alpha 15 Results\TMotor_Relaxed_J0.1109_0.0005_440.mat')

deg_per_ts = valRPM.*(pi/30).*(180/pi).*valDELTIME;
vecPOS_R = [0:(length(CT_U) - 1)]'.*deg_per_ts;
vecPOS = mod(vecPOS_R,360);

CT_ext = repmat(CT_U((cutoff*npts + 1):end),n,1);

fs = 1/valDELTIME; % sampling rate
rps = valRPM/60; % rotation rate
T = length(CT_ext).*valDELTIME; % total time
N = fs*T; % total samples
t = linspace(0, T, N);

FFTParam.FFTLength = fs; % notice how your BR peaks change when you adjust fft length here, e.g. fs*10
FFTParam.WinType = 'rectwin'; % play with flattopwin or rectwin to see difference windowing makes on spectra
FFTParam.SampRate = fs; 
FFTParam.DFlag = 'none'; 
FFTParam.PercOverlap = 50;

InputData = detrend(CT_ext);
[Pxy,COH, F, Navg] = CalculatePSDandCohs(InputData, FFTParam);

figure(2);
hold on
plot(F, 10*log10(Pxy), '-.b')
hold off

xlabel('Frequency (Hz)')
ylabel('Power/Frequency (dB/Hz)')
grid minor
box on
axis tight

title('Alpha 15, \mu = 0.1109')
xlim([0 1000])













% WH = [4.5 5];
% % fcnFIG2LATEX(hFig68, 'tmotor_fft_0.pdf', WH)





