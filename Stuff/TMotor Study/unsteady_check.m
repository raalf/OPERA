clc
clear

ds = 80

%% Getting tunnel data
load('G:\GIT\opera\Stuff\TMotor Study\Tunnel Testing\2019-11-06\06-Nov-2019 12.56.34_Scorpion_ASI_T-Motor 18in_RPM3000_Alpha15_9.8342.mat')
close all

%%

hFig99 = figure(99);
clf(99);

[vecPOS_TUNNEL_OG, idx] = sort(vecPOS_TUNNEL_OG, 'ascend');
CT_tunnel = CT_tunnel(idx);
plot(downsample(vecPOS_TUNNEL_OG, ds), downsample(CT_tunnel, ds), '--r')
% plot(vecPOS_TUNNEL_OG, CT_tunnel, '--r')
grid minor
box on
axis tight

%% Getting tunnel data
load('G:\GIT\opera\Stuff\TMotor Study\Tunnel Testing\2019-10-31\31-Oct-2019 12.41.35_Scorpion_ASI_T-Motor 18in_RPM3000_Alpha30_21.5736.mat', 'vecPOS_TUNNEL_OG', 'CT_tunnel')

%%
figure(99);
hold on
plot(downsample(vecPOS_TUNNEL_OG, ds), downsample(CT_tunnel, ds), '-b')
% plot(vecPOS_TUNNEL_OG, CT_tunnel, '-b')
hold off

legend('\mu = 0.1346', '\mu = 0.3004','Location','NorthWest')

xlabel('Azimuth Location, Degrees')
ylabel('Thrust Coefficient') 

fcnFIG2LATEX(hFig99, 'TMotor_tunnel_transient.pdf', [8 5])