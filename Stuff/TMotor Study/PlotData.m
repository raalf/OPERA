clear

load('Tunnel Testing\2020-01-30\30-Jan-2020 16.15.50_Scorpion_KDE_T-Motor 18in_RPM3000_Alpha0_20.7492.mat', 'CT_tunnel', 'valRPM', 'dataRate'); % Alpha 0, 0.2889

fs = dataRate; % sampling rate
rps = valRPM/60; % rotation rate
T = 4*60; % total time
N = fs*T; % total samples
t = linspace(0, T, N);

% Frequency Domain
FFTParam.FFTLength = fs; % notice how your BR peaks change when you adjust fft length here, e.g. fs*10
FFTParam.WinType = 'hanning'; % play with flattopwin or rectwin to see difference windowing makes on spectra

FFTParam.SampRate = fs; FFTParam.DFlag = 'none'; FFTParam.PercOverlap = 50;
InputData = CT_tunnel;
[Pxy,COH, F, Navg]=CalculatePSDandCohs(InputData, FFTParam);

hFig10 = figure(10);
clf(10);
plot(F, 10*log10(Pxy))
xlabel('Frequency, Hz'); 
ylabel('Autospectrum, dB')
box on



load('Tunnel Testing\2020-01-30\30-Jan-2020 16.40.19_Scorpion_KDE_T-Motor 18in_RPM3000_Alpha15_7.9674.mat', 'CT_tunnel', 'valRPM', 'dataRate'); % Alpha 15, 0.1109
InputData = CT_tunnel;
[Pxy,COH, F, Navg]=CalculatePSDandCohs(InputData, FFTParam);

yyaxis right
plot(F, 10*log10(Pxy))
ylabel('Autospectrum, dB')
grid minor
box on