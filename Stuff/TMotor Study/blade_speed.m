clc
clear

% load('Alpha 0 Results/TMotor_Relaxed_J0.3_0.0005.mat')
load('Alpha 15 Results/TMotor_Fixed_J0.1115_0.0005.mat')
load('Alpha 45 Results/TMotor_Fixed_J0.1109_0.0005.mat')

% idx1 = vecDVESURFACE == 1;
% idx2 = vecDVESURFACE == 2;
% 
% blade_vel(:,1) = sum(sqrt(sum(matUINF_ALL(idx1,:,:).^2,2)).*vecDVEAREA(idx1),1);
% blade_vel(:,2) = sum(sqrt(sum(matUINF_ALL(idx2,:,:).^2,2)).*vecDVEAREA(idx2),1);
% 
% 
% deg_per_ts = valRPM.*(pi/30).*(180/pi).*valDELTIME;
% tmp_offset = 0;
% vecPOS_R = [tmp_offset:(length(CT_U)-1 + tmp_offset)]'.*deg_per_ts + 90;
% 
% tmp = vecPOS_R./360;
% start = 2;
% 
% idx = tmp >= start & tmp < start+1;
% idx2 = tmp >= start+1 & tmp < start+2;
% CT_relaxed2 = CT(~isnan(CT));
% tmp2 = find(idx);
% offset = -vecPOS_R(tmp2(1));
% 
% hFig1 = figure(1);
% clf(1)
% hold on
% plot(vecPOS_R(idx) + offset, (CT_U(idx) + CT_U(idx2))./2, '--b');
% plot(vecPOS_R(idx) + offset, (CT(idx) + CT(idx2))./2, '-.r');
% hold off
% grid minor
% box on
% axis tight
% ylabel('Thrust Coefficient');
% 
% yyaxis right
% hold on
% % plot(vecPOS_R(idx) + offset, (blade_vel(idx,1) + blade_vel(idx2,1))./2, '-.k');
% % plot(vecPOS_R(idx) + offset, (blade_vel(idx,2) + blade_vel(idx2,2))./2, '-m');
% % plot(vecPOS_R(idx) + offset, (blade_vel(idx,1) + blade_vel(idx2,1) + blade_vel(idx,2) + blade_vel(idx2,2))./4, '--g');
% plot(vecPOS_R(idx) + offset, sum(matINTCIRC(idx,:),2), '-m');
% hold off
% xlabel('Azimuth Location (Degrees)');




