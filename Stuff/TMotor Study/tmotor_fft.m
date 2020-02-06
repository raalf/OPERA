clc
clear

load('Tunnel Testing\2020-01-30\30-Jan-2020 16.15.50_Scorpion_KDE_T-Motor 18in_RPM3000_Alpha0_20.7492.mat', ... % Alpha 0, 0.2889
    'Angle', 'lbf_N', 'FT', 'rho', 'valDIAM', 'valRPM', 'vecPOS_TUNNEL_OG', 'dataRate');

close all

% Tunnel
CT_tunnel = lbf_N.*FT(:,3);
CT_tunnel = CT_tunnel./(rho.*(pi.*((valDIAM/2).^2)).*(((valDIAM/2).*(valRPM.*(pi/30))).^2));

[pxx, f] = periodogram(CT_tunnel, [], [], dataRate);

plot(f,10*log10(pxx))

xlabel('Frequency (Hz)')
ylabel('Power/Frequency (dB/Hz)')

grid minor
box on
axis tight

% WH = [4.5 5];
% % fcnFIG2LATEX(hFig68, 'tmotor_fft_0.pdf', WH)

load('Tunnel Testing\2020-01-30\30-Jan-2020 16.40.19_Scorpion_KDE_T-Motor 18in_RPM3000_Alpha15_7.9674.mat', ... % Alpha 15, 0.1109
    'Angle', 'lbf_N', 'FT', 'rho', 'valDIAM', 'valRPM', 'vecPOS_TUNNEL_OG', 'dataRate');

CT_tunnel = lbf_N.*FT(:,3);
CT_tunnel = CT_tunnel./(rho.*(pi.*((valDIAM/2).^2)).*(((valDIAM/2).*(valRPM.*(pi/30))).^2));

[pxx, f] = periodogram(CT_tunnel, [], [], dataRate);
 hold on
plot(f,10*log10(pxx))
hold off
% xlabel('Frequency (Hz)')
% ylabel('Power/Frequency (dB/Hz)')

