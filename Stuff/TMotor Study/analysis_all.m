clc
clear

%% 0

%% 0.1080
load('Tunnel Testing\2020-02-11\11-Feb-2020 17.28.46_Scorpion_KDE_T-Motor 18in_RPM3000_Alpha0_7.7559.mat', ... % Alpha 0, 0.1080
    'Angle', 'lbf_N', 'FT', 'rho', 'valDIAM', 'valRPM', 'vecPOS_TUNNEL_OG', 'dataRate');
CT_tunnel_raw = lbf_N.*FT(:,3);
CT_tunnel_raw = CT_tunnel_raw./(rho.*(pi.*((valDIAM/2).^2)).*(((valDIAM/2).*(valRPM.*(pi/30))).^2));

binAng = linspace(0, 360, 30);
for i = 1:length(binAng)
    rng = 0.25;
    idx = vecPOS_TUNNEL_OG >= binAng(i) - rng & vecPOS_TUNNEL_OG <= binAng(i) + rng;
    
    binAvg(i) = mean(CT_tunnel_raw(idx));
    binMax(i) = max(CT_tunnel_raw(idx));
    binMin(i) = min(CT_tunnel_raw(idx));
end
CT(:,1) = binAng';
CT(:,2) = binAvg';

Fp = 1/dataRate;
Fst = 2000/dataRate;
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
CTF(:,1) = binAng';
CTF(:,2) = binAvg';

%% 0.1441
load('Tunnel Testing\2020-02-11\11-Feb-2020 17.24.34_Scorpion_KDE_T-Motor 18in_RPM3000_Alpha0_10.3457.mat', ... % Alpha 0, 0.1441
    'Angle', 'lbf_N', 'FT', 'rho', 'valDIAM', 'valRPM', 'vecPOS_TUNNEL_OG', 'dataRate');
CT_tunnel_raw = lbf_N.*FT(:,3);
CT_tunnel_raw = CT_tunnel_raw./(rho.*(pi.*((valDIAM/2).^2)).*(((valDIAM/2).*(valRPM.*(pi/30))).^2));

binAng = linspace(0, 360, 30);
for i = 1:length(binAng)
    rng = 0.25;
    idx = vecPOS_TUNNEL_OG >= binAng(i) - rng & vecPOS_TUNNEL_OG <= binAng(i) + rng;
    
    binAvg(i) = mean(CT_tunnel_raw(idx));
    binMax(i) = max(CT_tunnel_raw(idx));
    binMin(i) = min(CT_tunnel_raw(idx));
end
CT(:,3) = binAvg';

Fp = 1/dataRate;
Fst = 2000/dataRate;
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
CTF(:,3) = binAvg';

%% 0.2043
load('Tunnel Testing\2020-02-11\11-Feb-2020 17.18.43_Scorpion_KDE_T-Motor 18in_RPM3000_Alpha0_14.6752.mat', ... % Alpha 0, 0.2043
    'Angle', 'lbf_N', 'FT', 'rho', 'valDIAM', 'valRPM', 'vecPOS_TUNNEL_OG', 'dataRate');
CT_tunnel_raw = lbf_N.*FT(:,3);
CT_tunnel_raw = CT_tunnel_raw./(rho.*(pi.*((valDIAM/2).^2)).*(((valDIAM/2).*(valRPM.*(pi/30))).^2));

binAng = linspace(0, 360, 30);
for i = 1:length(binAng)
    rng = 0.25;
    idx = vecPOS_TUNNEL_OG >= binAng(i) - rng & vecPOS_TUNNEL_OG <= binAng(i) + rng;
    
    binAvg(i) = mean(CT_tunnel_raw(idx));
    binMax(i) = max(CT_tunnel_raw(idx));
    binMin(i) = min(CT_tunnel_raw(idx));
end
CT(:,4) = binAvg';

Fp = 1/dataRate;
Fst = 2000/dataRate;
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
CTF(:,4) = binAvg';

%% 0.2889
load('Tunnel Testing\2020-01-30\30-Jan-2020 16.15.50_Scorpion_KDE_T-Motor 18in_RPM3000_Alpha0_20.7492.mat', ... % Alpha 0, 0.2889
    'Angle', 'lbf_N', 'FT', 'rho', 'valDIAM', 'valRPM', 'vecPOS_TUNNEL_OG', 'dataRate');
CT_tunnel_raw = lbf_N.*FT(:,3);
CT_tunnel_raw = CT_tunnel_raw./(rho.*(pi.*((valDIAM/2).^2)).*(((valDIAM/2).*(valRPM.*(pi/30))).^2));

binAng = linspace(0, 360, 30);
for i = 1:length(binAng)
    rng = 0.25;
    idx = vecPOS_TUNNEL_OG >= binAng(i) - rng & vecPOS_TUNNEL_OG <= binAng(i) + rng;
    
    binAvg(i) = mean(CT_tunnel_raw(idx));
    binMax(i) = max(CT_tunnel_raw(idx));
    binMin(i) = min(CT_tunnel_raw(idx));
end
CT(:,5) = binAvg';

Fp = 1/dataRate;
Fst = 2000/dataRate;
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
CTF(:,5) = binAvg';

%%

% yl = [-0.01 0.035];
% xl = [0 360];

hFig5 = figure(5);
clf(5);

hold on
plot(CT(:,1), CT(:,2), '-k'); 
plot(CT(:,1), CT(:,3), '-.r'); 
plot(CT(:,1), CT(:,4), '--m'); 
plot(CT(:,1), CT(:,5), ':b'); 
hold off

box on
grid minor
axis tight
legend('\mu_{\infty} = 0.1080','\mu_{\infty} = 0.1441','\mu_{\infty} = 0.2043','\mu_{\infty} = 0.2889','Location','North')

xlabel('Angle (Degrees)');
ylabel('Thrust Coefficient');
% 
% xlim(xl);
% ylim(yl);

WH = [4.5 5];
fcnFIG2LATEX(hFig5, 'tmotor_time_02.pdf', WH)

%%
hFig6 = figure(6);
clf(6);

hold on
plot(CTF(:,1), CTF(:,2), '-k'); 
plot(CTF(:,1), CTF(:,3), '-.r'); 
plot(CTF(:,1), CTF(:,4), '--m'); 
plot(CTF(:,1), CTF(:,5), ':b'); 
hold off

box on
grid minor
axis tight
legend('\mu_{\infty} = 0.1080','\mu_{\infty} = 0.1441','\mu_{\infty} = 0.2043','\mu_{\infty} = 0.2889','Location','North')

xlabel('Angle (Degrees)');
ylabel('Thrust Coefficient');

% xlim(xl);
% ylim(yl);

WH = [4.5 5];
fcnFIG2LATEX(hFig6, 'tmotor_time_02F.pdf', WH)

