clc
clear

ds = 80

%% Getting tunnel data
load('G:\GIT\opera\Stuff\TMotor Study\Tunnel Testing\2019-11-06\06-Nov-2019 12.56.34_Scorpion_ASI_T-Motor 18in_RPM3000_Alpha15_9.8342.mat')
close all

hFig98 = figure(98);
clf(98);

binAng = linspace(0, 360, 30);
for i = 1:length(binAng)
    rng = 0.25;
    idx = vecPOS_TUNNEL_OG >= binAng(i) - rng & vecPOS_TUNNEL_OG <= binAng(i) + rng;
    
    binAvg(i) = mean(CT_tunnel(idx));
    binMax(i) = max(CT_tunnel(idx));
    binMin(i) = min(CT_tunnel(idx));
end
errorbar(binAng, binAvg, binMin - binAvg, binMax - binAvg, '--ks');
grid minor
box on
axis tight

%% Getting tunnel data
load('G:\GIT\opera\Stuff\TMotor Study\Tunnel Testing\2019-10-31\31-Oct-2019 12.41.35_Scorpion_ASI_T-Motor 18in_RPM3000_Alpha30_21.5736.mat') % 0.3004
close(1)
figure(98)
binAng = linspace(0, 360, 30);
for i = 1:length(binAng)
    rng = 0.25;
    idx = vecPOS_TUNNEL_OG >= binAng(i) - rng & vecPOS_TUNNEL_OG <= binAng(i) + rng;
    
    binAvg(i) = mean(CT_tunnel(idx));
    binMax(i) = max(CT_tunnel(idx));
    binMin(i) = min(CT_tunnel(idx));
end
hold on
errorbar(binAng, binAvg, binMin - binAvg, binMax - binAvg, '-.bo');
hold off

legend('\mu = 0.1346', '\mu = 0.3004','Location','NorthWest')

xlabel('Azimuth Location, Degrees')
ylabel('Thrust Coefficient') 

fcnFIG2LATEX(hFig98, 'TMotor_tunnel_transient.pdf', [8 5])