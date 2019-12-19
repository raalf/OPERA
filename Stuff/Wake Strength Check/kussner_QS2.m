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

hFig39 = figure(39);
clf(39);
subplot(1,2,1);
plot(((s(1:40))), c_l(1:40), '-k');
box on
axis tight
grid minor

hold on
% disp('THIS IS WITH APPARENT MASS!!!!!!!!!!')
load('QS2\m2_dt1.mat', 'CL', 'valDELTIME', 'valMAXTIME', 'valSPAN', 'valAREA');
AR = (valSPAN.^2)./valAREA;
CLOG2D = CL.*((AR + 2)/AR);
CL_QS2{:,1} = CLOG2D;
s_t = [1:valMAXTIME].*valDELTIME.*2;
s_QS2{:,1} = s_t;
plot(s_t((1/valDELTIME):end) - 2, CLOG2D((1/valDELTIME):end), 'ok');
M(1) = 2;
dt(1) = valDELTIME;
ratio(1) = valDELTIME/(1/M(1));

load('QS2\m5_dt1.mat', 'CL', 'valDELTIME', 'valMAXTIME', 'valSPAN', 'valAREA');
AR = (valSPAN.^2)./valAREA;
CLOG2D = CL.*((AR + 2)/AR);
CL_QS2{:,2} = CLOG2D;
s_t = [1:valMAXTIME].*valDELTIME.*2;
s_QS2{:,2} = s_t;
plot(s_t((1/valDELTIME):end) - 2, CLOG2D((1/valDELTIME):end), 'b^');
M(2) = 5;
dt(2) = valDELTIME;
ratio(2) = valDELTIME/(1/M(2));

load('QS2\m8_dt1.mat', 'CL', 'valDELTIME', 'valMAXTIME', 'valSPAN', 'valAREA');
AR = (valSPAN.^2)./valAREA;
CLOG2D = CL.*((AR + 2)/AR);
CL_QS2{:,3} = CLOG2D;
s_t = [1:valMAXTIME].*valDELTIME.*2;
s_QS2{:,3} = s_t;
plot(s_t((1/valDELTIME):end) - 2, CLOG2D((1/valDELTIME):end), 'rs');
M(3) = 8;
dt(3) = valDELTIME;
ratio(3) = valDELTIME/(1/M(3));
hold off

xlabel('Distanced Travelled in Semi-Chords');
ylabel('Lift Coefficient');

legend('Kussner Function', ...
    ['M = ', num2str(M(1)), ', \Deltax_w/\Deltax_c = ', num2str(ratio(1))], ...
    ['M = ', num2str(M(2)), ', \Deltax_w/\Deltax_c = ', num2str(ratio(2))], ...
    ['M = ', num2str(M(3)), ', \Deltax_w/\Deltax_c = ', num2str(ratio(3))], 'Location','SouthEast')

xl = xlim;
yl = ylim;

title('Coarse Timesteps (\Delta_T = 1)');

%%
subplot(1,2,2);
plot(((s(1:40))), c_l(1:40), '-k');
box on
axis tight
grid minor

hold on
load('QS2\m2_dt0.1.mat', 'CL', 'valDELTIME', 'valMAXTIME', 'valSPAN', 'valAREA');
AR = (valSPAN.^2)./valAREA;
CLOG2D = CL.*((AR + 2)/AR);
CL_QS2{:,4} = CLOG2D;
s_t = [1:valMAXTIME].*valDELTIME.*2;
s_QS2{:,4} = s_t;
plot(s_t((1/valDELTIME):end) - 2, CLOG2D((1/valDELTIME):end), 'ok');
M(4) = 2;
dt(4) = valDELTIME;
ratio(4) = valDELTIME/(1/M(4));

load('QS2\m5_dt0.1.mat', 'CL', 'valDELTIME', 'valMAXTIME', 'valSPAN', 'valAREA');
AR = (valSPAN.^2)./valAREA;
CLOG2D = CL.*((AR + 2)/AR);
CL_QS2{:,5} = CLOG2D;
s_t = [1:valMAXTIME].*valDELTIME.*2;
s_QS2{:,5} = s_t;
plot(s_t((1/valDELTIME):end) - 2, CLOG2D((1/valDELTIME):end), 'b^');
M(5) = 5;
dt(5) = valDELTIME;
ratio(5) = valDELTIME/(1/M(5));

load('QS2\m8_dt0.1.mat', 'CL', 'valDELTIME', 'valMAXTIME', 'valSPAN', 'valAREA');
AR = (valSPAN.^2)./valAREA;
CLOG2D = CL.*((AR + 2)/AR);
CL_QS2{:,6} = CLOG2D;
s_t = [1:valMAXTIME].*valDELTIME.*2;
s_QS2{:,6} = s_t;
plot(s_t((1/valDELTIME):end) - 2, CLOG2D((1/valDELTIME):end), 'rs');
M(6) = 8;
dt(6) = valDELTIME;
ratio(6) = valDELTIME/(1/M(6));

hold off
xlabel('Distanced Travelled in Semi-Chords');
ylabel('Lift Coefficient');

legend('Kussner Function', ...
    ['M = ', num2str(M(4)), ', \Deltax_w/\Deltax_c = ', num2str(ratio(4))], ...
    ['M = ', num2str(M(5)), ', \Deltax_w/\Deltax_c = ', num2str(ratio(5))], ...
    ['M = ', num2str(M(6)), ', \Deltax_w/\Deltax_c = ', num2str(ratio(6))], 'Location','SouthEast')

xlim(xl);
ylim(yl);

title('Fine Timesteps (\Delta_T = 0.1)');

sgtitle('Quasi-Steady II Kussner Response')
