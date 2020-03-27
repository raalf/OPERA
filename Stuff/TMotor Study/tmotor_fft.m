clc
clear

%% Tunnel
load('Tunnel Testing\2020-01-30\30-Jan-2020 16.15.50_Scorpion_KDE_T-Motor 18in_RPM3000_Alpha0_20.7492.mat', ... % Alpha 0, 0.2889
'Angle', 'lbf_N', 'FT', 'rho', 'valDIAM', 'valRPM', 'vecPOS_TUNNEL_OG', 'dataRate');

fs = dataRate; % sampling rate

FFTParam.FFTLength = fs; % notice how your BR peaks change when you adjust fft length here, e.g. fs*10
FFTParam.WinType = 'flattopwin'; % play with flattopwin or rectwin to see difference windowing makes on spectra
FFTParam.SampRate = fs; 
FFTParam.DFlag = 'none'; 
FFTParam.PercOverlap = 50;

% Unfiltered
CT_tunnel = lbf_N.*FT(:,3);
CT_tunnel = CT_tunnel./(rho.*(pi.*((valDIAM/2).^2)).*(((valDIAM/2).*(valRPM.*(pi/30))).^2));
InputData = detrend(CT_tunnel);
[Pxy,COH, F, Navg]=CalculatePSDandCohs(InputData, FFTParam);

hFig5 = figure(5);
clf(5);
plot(F, 10*log10(Pxy), '-k')

Fp = 1/dataRate;
Fst = 2000/(dataRate/2);
d = fdesign.lowpass('N,Fc,Ap,Ast',1000,Fst,1,50);

% Hd = design(d);
% CT_tunnel2 = filtfilt(Hd.Numerator,1,detrend(CT_tunnel)) + mean(CT_tunnel);
% InputData = detrend(CT_tunnel2);
% [Pxy,COH, F, Navg]=CalculatePSDandCohs(InputData, FFTParam);
% hold on
% plot(F, 10*log10(Pxy), '--b')
% hold off
% 
% % Cutoff line
% yl = ylim;
% hold on
% plot([2000 2000],ylim,'-r','linewidth',1)
% hold off

xlabel('Frequency, Hz')
ylabel('Force Autospectrum, dB re 1 N^2/hz')
% grid minor
% box on

WH = [4.5*2 5];
fcnFIG2LATEX(hFig5, 'correct.pdf', WH)

%%
% Alpha 15
load('G:\GIT\opera\Stuff\TMotor Study\Alpha 15 Results\New\TMotor_Fixed_J0.2113_0.00025_newint.mat')
CT_U = CT;
% load('G:\GIT\opera\Stuff\TMotor Study\Alpha 15 Results\New\TMotor_Relaxed_J0.2113_0.00025.mat')

deg_per_ts = valRPM.*(pi/30).*(180/pi).*valDELTIME;
vecPOS_R = [0:(length(CT_U) - 1)]'.*deg_per_ts;
vecPOS = mod(vecPOS_R,360);

npts = 80;
cutoff = 1;
n = ceil((4*60)*((valRPM/60)./(3 - cutoff))); % 4 minutes?
CT_ext = repmat(CT_U((cutoff*npts + 1):end),n,1);

% hFig2 = figure(2)
% clf(2);
% plot(CT_ext)

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

hold on
plot(F, 10*log10(Pxy), '-.m')
hold off

grid minor
box on
axis tight
title('Alpha 0, \mu = 0.3')

xlim([0 1000])

legend('Experimental','Fixed Wake (80 Az/rev)','Location','NorthEast')
% % Alpha 15
% load('G:\GIT\opera\Stuff\TMotor Study\Alpha 15 Results\TMotor_Relaxed_J0.1109_0.0005_440.mat')
% 
% deg_per_ts = valRPM.*(pi/30).*(180/pi).*valDELTIME;
% vecPOS_R = [0:(length(CT_U) - 1)]'.*deg_per_ts;
% vecPOS = mod(vecPOS_R,360);
% 
% CT_ext = repmat(CT_U((cutoff*npts + 1):end),n,1);
% 
% fs = 1/valDELTIME; % sampling rate
% rps = valRPM/60; % rotation rate
% T = length(CT_ext).*valDELTIME; % total time
% N = fs*T; % total samples
% t = linspace(0, T, N);
% 
% FFTParam.FFTLength = fs; % notice how your BR peaks change when you adjust fft length here, e.g. fs*10
% FFTParam.WinType = 'rectwin'; % play with flattopwin or rectwin to see difference windowing makes on spectra
% FFTParam.SampRate = fs; 
% FFTParam.DFlag = 'none'; 
% FFTParam.PercOverlap = 50;
% 
% InputData = detrend(CT_ext);
% [Pxy,COH, F, Navg] = CalculatePSDandCohs(InputData, FFTParam);
% 
% figure(2);
% hold on
% plot(F, 10*log10(Pxy), '-.b')
% hold off
% 
% xlabel('Frequency (Hz)')
% ylabel('Power/Frequency (dB/Hz)')
% grid minor
% box on
% axis tight
% 
% title('Alpha 15, \mu = 0.1109')
% xlim([0 1000])




% % title('Alpha 15, \mu = 0.2113')
% % xlim([0 1000])
% 
% % legend('Experimental','DDE Method','Location','NorthEast');
% WH = [4.5*2 5];
% fcnFIG2LATEX(hFig2, 'tmotor_fft_15.pdf', WH)
% 
% blade_rate_table_15 = [100 -55.9845 -56.2687; 200 -56.8202 -73.7387; 300 -63.1627 -76.6803; 400 -76.9571 -74.3934];
% blade_rate_table_0 = [100 -47.5829 -51.7143; 200 -48.2955 -59.3588; 300 -54.6626 -67.4689; 400 -71.8582 -66.5435];
% 
% 
% 
% 
% 
% 
% 
% 
WH = [4.5*2 5];
fcnFIG2LATEX(hFig5, 'tmotor_autospectrum_15_0.2113.pdf', WH)





