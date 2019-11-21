clc
clear

%%
% Kussner function, numerical integration
% C_k = @(k)(besselh(1,2,k))./((besselh(1,2,k)) + 1i.*(besselh(0,2,k)));
% tmp = @(k) C_k(k).*(besselj(0,k) - 1i.*besselj(1,k)) + 1i.*besselj(1,k);
% F_k = @(k)real(tmp(k));
% G_k = @(k)imag(tmp(k));
% 
% 
% w0 = 0.055;
% vinf = 1;
% s = linspace(0,20,100)';
% 
% for i = 1:length(s)
%     tmp2 = @(k) ((F_k(k).*cos(k) + G_k(k).*sin(k))./k).*sin(k.*s(i));
%     phi_s = (2./pi).*integral(tmp2,0,1e4);
%     c_l(i) = (2.*pi).*(w0./vinf).*phi_s;
% end

load('Kussner.mat')

hFig22 = figure(22);
clf(22);
plot(s, c_l, '-k');
box on
axis tight
grid minor

hold on
load('Kussner_m2_0.5_fixed.mat', 'CL2D', 'valDELTIME', 'valMAXTIME', 'valGUSTSTART');
s_t = ([1:valMAXTIME-1] - valGUSTSTART).*valDELTIME.*2;
m2_err = abs(CL2D(2:end) - interp1(s,c_l,s_t)')./interp1(s,c_l,s_t)';

% figure(23);
% clf(23);
% plot(s_t, err, '-sb');
% box on
% axis tight
% grid minor

figure(22);
plot(s_t, CL2D(2:end), '-sb');

load('Kussner_m8_0.25_fixed.mat', 'CL2D', 'valDELTIME', 'valMAXTIME', 'valGUSTSTART');
s_t = ([1:valMAXTIME-1] - valGUSTSTART).*valDELTIME.*2;
m4_err = abs(CL2D(2:end) - interp1(s,c_l,s_t)')./interp1(s,c_l,s_t)';

% figure(23);
% hold on
% plot(s_t, err, '-.^r');
% hold off
% box on
% axis tight
% grid minor

figure(22);
plot(s_t, CL2D(2:end), '-.^r');

load('Kussner_m8_0.125_fixed.mat', 'CL2D', 'valDELTIME', 'valMAXTIME');
s_t = ([1:valMAXTIME-1] - valGUSTSTART).*valDELTIME.*2;
m8_err = abs(CL2D(2:end) - interp1(s,c_l,s_t)')./interp1(s,c_l,s_t)';

% figure(23);
% hold on
% plot(s_t, err, '--m');
% hold off
% box on
% axis tight
% grid minor

figure(22);
plot(s_t, CL2D(2:end), '--m');

hold off

legend('Kussner Function', 'M = 2', 'M = 4', 'M = 8', 'Location', 'SouthEast')
xlabel('Distance Travelled by Gust in Semi-Chords')
ylabel('Two-Dimensional Lift Coefficient')

% fcnFIG2LATEX(hFig22, 'Kussner.pdf', [8 5])


%%

T = table([2 4 8]', [1 1 1]', [nanmean(m2_err(~isinf(m2_err))); nanmean(m4_err(~isinf(m4_err))); nanmean(m8_err(~isinf(m8_err)))].*100, 'VariableNames', {'M', 'Distance_Chord_Ratio', 'Average_Percent_Error'})







