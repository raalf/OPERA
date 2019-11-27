clc
clear

%%
% % Kussner function, numerical integration
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

hFig25 = figure(25);
clf(25);
subplot(1,2,1);
plot(((s(1:40))), c_l(1:40), '-k');
box on
axis tight
grid minor

hold on
% disp('THIS IS WITH APPARENT MASS!!!!!!!!!!')
load('QS2\m2_dt1.mat');
[~, CL_U, ~] = fcnDGAMMADT(1, valDELTIME, strATYPE, matINTCIRC, valDENSITY, valRPM, valDIAM, valAREA, valUINF, matLIFTFREE, matLIFTIND, matDRAGIND, matSIDEFREE, matSIDEIND, matDVELIFT_DIR, matDVEDRAG_DIR, matDVESIDE_DIR);
s_t = [1:valMAXTIME].*valDELTIME.*2;
valAR = (valSPAN.^2)./valAREA;
CLOG2D = CL_U.*((valAR + 2)/valAR);
idx = (1/valDELTIME):length(s_t);
m2c_err = abs(CLOG2D(idx) - interp1(s,c_l,s_t(idx)-2)')./interp1(s,c_l,s_t(idx)-2)';
m2c_err = [s_t(idx)' - 2 m2c_err];
plot(s_t(idx) - 2, CLOG2D(idx), 'ok');
M(1) = 2;
ratio(1) = valDELTIME/(1/M(1));

load('QS2\m5_dt1.mat');
[~, CL_U, ~] = fcnDGAMMADT(1, valDELTIME, strATYPE, matINTCIRC, valDENSITY, valRPM, valDIAM, valAREA, valUINF, matLIFTFREE, matLIFTIND, matDRAGIND, matSIDEFREE, matSIDEIND, matDVELIFT_DIR, matDVEDRAG_DIR, matDVESIDE_DIR);
s_t = [1:valMAXTIME].*valDELTIME.*2;
valAR = (valSPAN.^2)./valAREA;
CLOG2D = CL_U.*((valAR + 2)/valAR);
idx = (1/valDELTIME):length(s_t);
m5c_err = abs(CLOG2D(idx) - interp1(s,c_l,s_t(idx)-2)')./interp1(s,c_l,s_t(idx)-2)';
m5c_err = [s_t(idx)' - 2 m5c_err];
plot(s_t(idx) - 2, CLOG2D(idx), 'b^');
M(2) = 5;
ratio(2) = valDELTIME/(1/M(2));

load('QS2\m8_dt1.mat');
[~, CL_U, ~] = fcnDGAMMADT(1, valDELTIME, strATYPE, matINTCIRC, valDENSITY, valRPM, valDIAM, valAREA, valUINF, matLIFTFREE, matLIFTIND, matDRAGIND, matSIDEFREE, matSIDEIND, matDVELIFT_DIR, matDVEDRAG_DIR, matDVESIDE_DIR);
s_t = [1:valMAXTIME].*valDELTIME.*2;
valAR = (valSPAN.^2)./valAREA;
CLOG2D = CL_U.*((valAR + 2)/valAR);
idx = (1/valDELTIME):length(s_t);
m8c_err = abs(CLOG2D(idx) - interp1(s,c_l,s_t(idx)-2)')./interp1(s,c_l,s_t(idx)-2)';
m8c_err = [s_t(idx)' - 2 m8c_err];
plot(s_t(idx) - 2, CLOG2D(idx), 'rs');
M(3) = 8;
ratio(3) = valDELTIME/(1/M(3));


load('QS2\m10_dt1.mat');
skip = 1;
[~, CL_U, ~] = fcnDGAMMADT(skip, valDELTIME, strATYPE, matINTCIRC, valDENSITY, valRPM, valDIAM, valAREA, valUINF, matLIFTFREE, matLIFTIND, matDRAGIND, matSIDEFREE, matSIDEIND, matDVELIFT_DIR, matDVEDRAG_DIR, matDVESIDE_DIR);
s_t = [1:skip:valMAXTIME].*valDELTIME.*2;
valAR = (valSPAN.^2)./valAREA;
CLOG2D = CL_U.*((valAR + 2)/valAR);
idx = ceil(1/(valDELTIME.*skip)):length(s_t);
m10c_err = abs(CLOG2D(idx) - interp1(s,c_l,s_t(idx)-2)')./interp1(s,c_l,s_t(idx)-2)';
m10c_err = [s_t(idx)' - 2 m10c_err];
plot(s_t(idx) - 2, CLOG2D(idx), 'md');
M(8) = 10;
ratio(8) = (valDELTIME*skip)/(1/M(8));

hold off

xlabel('Distanced Travelled in Semi-Chords');
ylabel('Lift Coefficient');

legend('Kussner Function', ...
    ['M = ', num2str(M(1)), ', \Deltax_w/\Deltax_c = ', num2str(ratio(1))], ...
    ['M = ', num2str(M(2)), ', \Deltax_w/\Deltax_c = ', num2str(ratio(2))], ...
    ['M = ', num2str(M(3)), ', \Deltax_w/\Deltax_c = ', num2str(ratio(3))], ...
    ['M = ', num2str(M(8)), ', \Deltax_w/\Deltax_c = ', num2str(ratio(8))], 'Location','SouthEast')
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
load('QS2\m2_dt0.1.mat');
skip = 4;
[~, CL_U, ~] = fcnDGAMMADT(skip, valDELTIME, strATYPE, matINTCIRC, valDENSITY, valRPM, valDIAM, valAREA, valUINF, matLIFTFREE, matLIFTIND, matDRAGIND, matSIDEFREE, matSIDEIND, matDVELIFT_DIR, matDVEDRAG_DIR, matDVESIDE_DIR);
s_t = [1:skip:valMAXTIME].*valDELTIME.*2;
valAR = (valSPAN.^2)./valAREA;
CLOG2D = CL_U.*((valAR + 2)/valAR);
idx = ceil(1/(valDELTIME.*skip)):length(s_t);
m2_err = abs(CLOG2D(idx) - interp1(s,c_l,s_t(idx)-2)')./interp1(s,c_l,s_t(idx)-2)';
m2_err = [s_t(idx)' - 2 m2_err];
plot(s_t(idx) - 2, CLOG2D(idx), 'ok');
M(4) = 2;
ratio(4) = (valDELTIME*skip)/(1/M(4));

load('QS2\m5_dt0.1.mat');
skip = 4;
[~, CL_U, ~] = fcnDGAMMADT(skip, valDELTIME, strATYPE, matINTCIRC, valDENSITY, valRPM, valDIAM, valAREA, valUINF, matLIFTFREE, matLIFTIND, matDRAGIND, matSIDEFREE, matSIDEIND, matDVELIFT_DIR, matDVEDRAG_DIR, matDVESIDE_DIR);
s_t = [1:skip:valMAXTIME].*valDELTIME.*2;
valAR = (valSPAN.^2)./valAREA;
CLOG2D = CL_U.*((valAR + 2)/valAR);
idx = ceil(1/(valDELTIME.*skip)):length(s_t);
m5_err = abs(CLOG2D(idx) - interp1(s,c_l,s_t(idx)-2)')./interp1(s,c_l,s_t(idx)-2)';
m5_err = [s_t(idx)' - 2 m5_err];
plot(s_t(idx) - 2, CLOG2D(idx), 'b^');
M(5) = 5;
ratio(5) = (valDELTIME*skip)/(1/M(5));

load('QS2\m8_dt0.1.mat');
skip = 4;
[~, CL_U, ~] = fcnDGAMMADT(skip, valDELTIME, strATYPE, matINTCIRC, valDENSITY, valRPM, valDIAM, valAREA, valUINF, matLIFTFREE, matLIFTIND, matDRAGIND, matSIDEFREE, matSIDEIND, matDVELIFT_DIR, matDVEDRAG_DIR, matDVESIDE_DIR);
s_t = [1:skip:valMAXTIME].*valDELTIME.*2;
valAR = (valSPAN.^2)./valAREA;
CLOG2D = CL_U.*((valAR + 2)/valAR);
idx = ceil(1/(valDELTIME.*skip)):length(s_t);
m8_err = abs(CLOG2D(idx) - interp1(s,c_l,s_t(idx)-2)')./interp1(s,c_l,s_t(idx)-2)';
m8_err = [s_t(idx)' - 2 m8_err];
plot(s_t(idx) - 2, CLOG2D(idx), 'rs');
M(6) = 8;
ratio(6) = (valDELTIME*skip)/(1/M(6));

load('QS2\m10_dt0.1.mat');
skip = 4;
[~, CL_U, ~] = fcnDGAMMADT(skip, valDELTIME, strATYPE, matINTCIRC, valDENSITY, valRPM, valDIAM, valAREA, valUINF, matLIFTFREE, matLIFTIND, matDRAGIND, matSIDEFREE, matSIDEIND, matDVELIFT_DIR, matDVEDRAG_DIR, matDVESIDE_DIR);
s_t = [1:skip:valMAXTIME].*valDELTIME.*2;
valAR = (valSPAN.^2)./valAREA;
CLOG2D = CL_U.*((valAR + 2)/valAR);
idx = ceil(1/(valDELTIME.*skip)):length(s_t);
m10_err = abs(CLOG2D(idx) - interp1(s,c_l,s_t(idx)-2)')./interp1(s,c_l,s_t(idx)-2)';
m10_err = [s_t(idx)' - 2 m10_err];
plot(s_t(idx) - 2, CLOG2D(idx), 'md');
M(7) = 10;
ratio(7) = (valDELTIME*skip)/(1/M(7));

hold off
xlabel('Distanced Travelled in Semi-Chords');
ylabel('Lift Coefficient');

legend('Kussner Function', ...
    ['M = ', num2str(M(4)), ', \Deltax_w/\Deltax_c = ', num2str(ratio(4))], ...
    ['M = ', num2str(M(5)), ', \Deltax_w/\Deltax_c = ', num2str(ratio(5))], ...
    ['M = ', num2str(M(6)), ', \Deltax_w/\Deltax_c = ', num2str(ratio(6))], ...
    ['M = ', num2str(M(7)), ', \Deltax_w/\Deltax_c = ', num2str(ratio(7))], 'Location','SouthEast')

xlim(xl);
ylim(yl);

title('Fine Timesteps (\Delta_T = 0.1)');

sgtitle('Fully Unsteady Kussner Response')

%%

% hFig27 = figure(27);
% clf(27);
% 
% load('QS2\m2_dt0.1.mat');
% s_t = [1:valMAXTIME].*valDELTIME.*2;
% idx = ceil(1/valDELTIME):length(s_t);
% plot(s_t(idx) - 2, sum(matINTCIRC(idx,:),2), '--ok')
% 
% box on
% axis tight
% grid minor
% 
% xlabel('Distanced Travelled in Semi-Chords');
% ylabel('Integrated Circulation on Surface');

%%

hFig70 = figure(70);
clf(70);

subplot(1,2,2)
plot(m2_err(:,1), m2_err(:,2).*100, '-ok');
hold on
plot(m5_err(:,1), m5_err(:,2).*100, '-.b^');
plot(m8_err(:,1), m8_err(:,2).*100, '--rs');
plot(m10_err(:,1), m10_err(:,2).*100, '-md');
hold off
grid minor
box on
axis tight
xlabel('Distanced Travelled in Semi-Chords');
ylabel('Percent Error in Kussner Response');

legend(['M = ', num2str(M(4)), ', \Deltax_w/\Deltax_c = ', num2str(ratio(4))], ...
    ['M = ', num2str(M(5)), ', \Deltax_w/\Deltax_c = ', num2str(ratio(5))], ...
    ['M = ', num2str(M(6)), ', \Deltax_w/\Deltax_c = ', num2str(ratio(6))], ...
    ['M = ', num2str(M(7)), ', \Deltax_w/\Deltax_c = ', num2str(ratio(7))], 'Location','NorthEast')
xlim(xl)
y2 = ylim;

title('Fine Timesteps (\Delta_T = 0.1)');

subplot(1,2,1)
plot(m2c_err(:,1), m2c_err(:,2).*100, '-ok');
hold on
plot(m5c_err(:,1), m5c_err(:,2).*100, '-.b^');
plot(m8c_err(:,1), m8c_err(:,2).*100, '--rs');
plot(m10c_err(:,1), m10c_err(:,2).*100, '-md');
hold off
grid minor
box on
axis tight
xlabel('Distanced Travelled in Semi-Chords');
ylabel('Percent Error in Kussner Response');

legend(['M = ', num2str(M(4)), ', \Deltax_w/\Deltax_c = ', num2str(ratio(4))], ...
    ['M = ', num2str(M(5)), ', \Deltax_w/\Deltax_c = ', num2str(ratio(5))], ...
    ['M = ', num2str(M(6)), ', \Deltax_w/\Deltax_c = ', num2str(ratio(6))], ...
    ['M = ', num2str(M(7)), ', \Deltax_w/\Deltax_c = ', num2str(ratio(7))], 'Location','NorthEast')
xlim(xl)
ylim(y2)

title('Coarse Timesteps (\Delta_T = 1)');

sgtitle('Error in the Fully Unsteady Kussner Response')