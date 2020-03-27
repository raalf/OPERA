clc
clear

Fst_hz = 2000

%% 0

% 0.1080 ----------------------------------------------------------------
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
    binStd(i) = std(CT_tunnel_raw(idx));
end
% CV = max(abs(binStd./binAvg))
% er = [abs(binMax - binAvg)./binAvg abs(binMin - binAvg)./binAvg]'
CT(:,1) = binAng';
CT(:,2) = binAvg';
SG(:,2) = binStd';

Fst = Fst_hz/(dataRate/2);
d = fdesign.lowpass('N,Fc,Ap,Ast',1000,Fst,1,50);
% Fp = 1/dataRate;
% Fst = 2000/dataRate;
% d = fdesign.lowpass('N,Fp,Fst',3,Fp,Fst);
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

% 0.1441 ----------------------------------------------------------------
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
    binStd(i) = std(CT_tunnel_raw(idx));
end
% CV = max(abs(binStd./binAvg))
CT(:,3) = binAvg';
SG(:,3) = binStd';

Fst = Fst_hz/(dataRate/2);
d = fdesign.lowpass('N,Fc,Ap,Ast',1000,Fst,1,50);
% Fp = 1/dataRate;
% Fst = 2000/dataRate;
% d = fdesign.lowpass('N,Fp,Fst',3,Fp,Fst);
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

% 0.2043 ----------------------------------------------------------------
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
    binStd(i) = std(CT_tunnel_raw(idx));
end
% CV = max(abs(binStd./binAvg))
CT(:,4) = binAvg';
SG(:,4) = binStd';

Fst = Fst_hz/(dataRate/2);
d = fdesign.lowpass('N,Fc,Ap,Ast',1000,Fst,1,50);
% Fp = 1/dataRate;
% Fst = 2000/dataRate;
% d = fdesign.lowpass('N,Fp,Fst',3,Fp,Fst);
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

% 0.2889 ----------------------------------------------------------------
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
    binStd(i) = std(CT_tunnel_raw(idx));
end
% CV = max(abs(binStd./binAvg))
CT(:,5) = binAvg';
SG(:,5) = binStd';

Fst = Fst_hz/(dataRate/2);
d = fdesign.lowpass('N,Fc,Ap,Ast',1000,Fst,1,50);
% Fp = 1/dataRate;
% Fst = 2000/dataRate;
% d = fdesign.lowpass('N,Fp,Fst',3,Fp,Fst);
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

%  ----------------------------------------------------------------

% yl = [-0.01 0.035];
% xl = [0 360];

hFig5 = figure(5);
clf(5);

hold on
plot(CT(:,1), CT(:,2), '-k'); 
plot(CT(:,1), CT(:,3), '-.r'); 
plot(CT(:,1), CT(:,4), '--m'); 
% plot(CT(:,1), CT(:,5), ':b'); 
cs = 3;
% errorbar(CT(:,1), CT(:,2), SG(:,2), '-k', 'CapSize', cs); 
% errorbar(CT(:,1), CT(:,3), SG(:,3), '-.r', 'CapSize', cs); 
% errorbar(CT(:,1), CT(:,4), SG(:,4), '--m', 'CapSize', cs); 
errorbar(CT(:,1), CT(:,5), SG(:,5), ':b', 'CapSize', cs);
hold off

box on
grid minor
axis tight
legend('\mu_{\infty} = 0.1080','\mu_{\infty} = 0.1441','\mu_{\infty} = 0.2043','\mu_{\infty} = 0.2889','Location','North','FontSize',7)

xlabel('Azimuth Angle, Degrees');
ylabel('Thrust Coefficient');

xl = xlim;
yl = ylim;

WH = [4.5 5];
fcnFIG2LATEX(hFig5, 'tmotor_time_02.pdf', WH)

%  ----------------------------------------------------------------
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
legend('\mu_{\infty} = 0.1080','\mu_{\infty} = 0.1441','\mu_{\infty} = 0.2043','\mu_{\infty} = 0.2889','Location','North','FontSize',7)

xlabel('Azimuth Angle, Degrees');
ylabel('Thrust Coefficient');

xlim(xl);
ylim(yl);

WH = [4.5 5];
fcnFIG2LATEX(hFig6, 'tmotor_time_02F.pdf', WH)

%% 15
mu15 = [0.0848 0.1109 0.1480 0.2113];
% 0.0848 ----------------------------------------------------------------
load('Tunnel Testing\2020-02-11\11-Feb-2020 16.41.03_Scorpion_KDE_T-Motor 18in_RPM3000_Alpha15_6.0928.mat', ... % Alpha 15, 0.0848
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
    binStd(i) = std(CT_tunnel_raw(idx));

end
CT15(:,1) = binAng';
CT15(:,2) = binAvg';
SG(:,2) = binStd';

Fst = Fst_hz/(dataRate/2);
d = fdesign.lowpass('N,Fc,Ap,Ast',1000,Fst,1,50);
% Fp = 1/dataRate;
% Fst = 2000/dataRate;
% d = fdesign.lowpass('N,Fp,Fst',3,Fp,Fst);
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
CTF15(:,1) = binAng';
CTF15(:,2) = binAvg';

% 0.1109 ----------------------------------------------------------------
load('Tunnel Testing\2020-01-30\30-Jan-2020 16.40.19_Scorpion_KDE_T-Motor 18in_RPM3000_Alpha15_7.9674.mat', ... % Alpha 15, 0.1109
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
    binStd(i) = std(CT_tunnel_raw(idx));
end
CT15(:,3) = binAvg';
SG(:,3) = binStd';

Fst = Fst_hz/(dataRate/2);
d = fdesign.lowpass('N,Fc,Ap,Ast',1000,Fst,1,50);
% Fp = 1/dataRate;
% Fst = 2000/dataRate;
% d = fdesign.lowpass('N,Fp,Fst',3,Fp,Fst);
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
CTF15(:,3) = binAvg';

% 0.1480 ----------------------------------------------------------------
load('Tunnel Testing\2020-02-11\11-Feb-2020 17.01.50_Scorpion_KDE_T-Motor 18in_RPM3000_Alpha15_10.6312.mat', ... % Alpha 15, 0.1480
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
    binStd(i) = std(CT_tunnel_raw(idx));
end
CT15(:,4) = binAvg';
SG(:,4) = binStd';

Fst = Fst_hz/(dataRate/2);
d = fdesign.lowpass('N,Fc,Ap,Ast',1000,Fst,1,50);
% Fp = 1/dataRate;
% Fst = 2000/dataRate;
% d = fdesign.lowpass('N,Fp,Fst',3,Fp,Fst);
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
CTF15(:,4) = binAvg';

% 0.2113 ----------------------------------------------------------------
load('Tunnel Testing\2020-02-11\11-Feb-2020 17.06.36_Scorpion_KDE_T-Motor 18in_RPM3000_Alpha15_15.175.mat', ... % Alpha 15, 0.2113
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
    binStd(i) = std(CT_tunnel_raw(idx));
end
CT15(:,5) = binAvg';
SG(:,5) = binStd';

Fst = Fst_hz/(dataRate/2);
d = fdesign.lowpass('N,Fc,Ap,Ast',1000,Fst,1,50);
% Fp = 1/dataRate;
% Fst = 2000/dataRate;
% d = fdesign.lowpass('N,Fp,Fst',3,Fp,Fst);
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
    binStd(i) = std(CT_tunnel_raw(idx));
end
CTF15(:,5) = binAvg';

%  ----------------------------------------------------------------

hFig5 = figure(5);
clf(5);

hold on
plot(CT15(:,1), CT15(:,2), '-k'); 
plot(CT15(:,1), CT15(:,3), '-.r'); 
plot(CT15(:,1), CT15(:,4), '--m'); 
% plot(CT15(:,1), CT15(:,5), ':b'); 
cs = 3;
% errorbar(CT15(:,1), CT15(:,2), SG(:,2), '-k', 'CapSize', cs); 
% errorbar(CT15(:,1), CT15(:,3), SG(:,3), '-.r', 'CapSize', cs); 
% errorbar(CT15(:,1), CT15(:,4), SG(:,4), '--m', 'CapSize', cs); 
errorbar(CT15(:,1), CT15(:,5), SG(:,5), ':b', 'CapSize', cs);
hold off

box on
grid minor
axis tight
% 0.0848 0.1109 0.1480 0.2113
legend('\mu_{\infty} = 0.0848','\mu_{\infty} = 0.1109','\mu_{\infty} = 0.1480','\mu_{\infty} = 0.2113','Location','North','FontSize',7)

xlabel('Azimuth Angle, Degrees');
ylabel('Thrust Coefficient');

xl = xlim;
yl = ylim;

WH = [4.5 5];
fcnFIG2LATEX(hFig5, 'tmotor_time_152.pdf', WH)

%  ----------------------------------------------------------------
hFig6 = figure(6);
clf(6);

hold on
plot(CTF15(:,1), CTF15(:,2), '-k'); 
plot(CTF15(:,1), CTF15(:,3), '-.r'); 
plot(CTF15(:,1), CTF15(:,4), '--m'); 
plot(CTF15(:,1), CTF15(:,5), ':b'); 
hold off

box on
grid minor
axis tight
legend('\mu_{\infty} = 0.0848','\mu_{\infty} = 0.1109','\mu_{\infty} = 0.1480','\mu_{\infty} = 0.2113','Location','North','FontSize',7)

xlabel('Azimuth Angle, Degrees');
ylabel('Thrust Coefficient');

xlim(xl);
ylim(yl);

WH = [4.5 5];
fcnFIG2LATEX(hFig6, 'tmotor_time_152F.pdf', WH)

%%

% 0
for i = 2:5
    [pk,lk] = findpeaks(CT(:,i),CT(:,1),'MinPeakDistance',140,'SortStr','descend');
    [pk2,lk2] = findpeaks(-CT(:,i),CT(:,1),'MinPeakDistance',140,'SortStr','descend');
    CTA(i-1) = (mean(pk(1:2)) - mean(-pk2(1:2)))./mean(CT(:,i))./2;
end

for i = 2:5
    [pk,lk] = findpeaks(CTF(:,i),CTF(:,1),'MinPeakDistance',140,'SortStr','descend');
    [pk2,lk2] = findpeaks(-CTF(:,i),CTF(:,1),'MinPeakDistance',140,'SortStr','descend');
    CTFA(i-1) = (mean(pk(1:2)) - mean(-pk2(1:2)))./mean(CTF(:,i))./2;
end

% 15
for i = 2:5
    [pk,lk] = findpeaks(CT15(:,i),CT15(:,1),'MinPeakDistance',140,'SortStr','descend');
    [pk2,lk2] = findpeaks(-CT15(:,i),CT15(:,1),'MinPeakDistance',140,'SortStr','descend');
    CTA15(i-1) = (mean(pk(1:2)) - mean(-pk2(1:2)))./mean(CT15(:,i))./2;
end

for i = 2:5
    [pk,lk] = findpeaks(CTF15(:,i),CTF15(:,1),'MinPeakDistance',140,'SortStr','descend');
    [pk2,lk2] = findpeaks(-CTF15(:,i),CTF15(:,1),'MinPeakDistance',140,'SortStr','descend');
    CTFA15(i-1) = (mean(pk(1:2)) - mean(-pk2(1:2)))./mean(CTF15(:,i))./2;
end

mu = [0.1080 0.1441 0.2043 0.2889];

hFig7 = figure(7);
clf(7);
plot(mu, CTA, 'sk');
hold on
plot(mu, CTFA, 'bo');
plot(mu15, CTA15, 'm^');
plot(mu15, CTFA15, 'rd');
hold off
% lsline


hold on
[x_fit,y_fit] = quickfit(mu',CTA');
plot(x_fit,y_fit,'-k')

[x_fit,y_fit] = quickfit(mu',CTFA');
plot(x_fit,y_fit,'-b')

[x_fit,y_fit] = quickfit(mu15',CTA15');
plot(x_fit,y_fit,'-m')

[x_fit,y_fit] = quickfit(mu15',CTFA15');
plot(x_fit,y_fit,'-r')
hold off

grid minor
box on
axis tight
legend(['\alpha_{tpp} = 0', char(176)], ['\alpha_{tpp} = 0', char(176), ' (Low-Pass Filtered)'], ['\alpha_{tpp} = 15', char(176)], ['\alpha_{tpp} = 15', char(176), ' (Low-Pass Filtered)'], 'Location', 'NorthWest')
ylabel('Normalized Thrust Oscillation Amplitude')
xlabel('Normalized Freestream Velocity, V/(\OmegaR)')

xlim([0 0.3])
ylim([0 max(ylim)])
WH = [4.5*2 5];
fcnFIG2LATEX(hFig7, 'tmotor_oscill.pdf', WH)



function [x_fit,y_fit] = quickfit(mu,ct)
    ft1 = fittype({'x'});
    p1 = fit(mu,ct,ft1);
    x_fit = linspace(0,0.3,2);
    y_fit = feval(p1, x_fit);
end























