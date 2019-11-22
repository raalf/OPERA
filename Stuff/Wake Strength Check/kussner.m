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
plot(((s+2)/2), c_l, '-k', 'linewidth',2);
box on
axis tight
grid minor

hold on
% disp('THIS IS WITH APPARENT MASS!!!!!!!!!!')
load('m5_dt1_wfix2.mat');
[~, CL_U, ~] = fcnDGAMMADT(1, valDELTIME, strATYPE, matINTCIRC, valDENSITY, valRPM, valDIAM, valAREA, valUINF, matLIFTFREE, matLIFTIND, matDRAGIND, matSIDEFREE, matSIDEIND, matDVELIFT_DIR, matDVEDRAG_DIR, matDVESIDE_DIR);
s_t = [1:valMAXTIME].*valDELTIME;
valAR = (valSPAN.^2)./valAREA;
CL_U2D = CL_U.*((valAR + 2)/valAR);
plot(s_t, CL_U2D, '-ok', 'LineWidth', 1);
delt_time(1) = valDELTIME;
ratio(1) = valDELTIME/(1/5);

% load('m5_dt0.5_wfix2.mat');
% [~, CL_U, ~] = fcnDGAMMADT(1, valDELTIME, strATYPE, matINTCIRC, valDENSITY, valRPM, valDIAM, valAREA, valUINF, matLIFTFREE, matLIFTIND, matDRAGIND, matSIDEFREE, matSIDEIND, matDVELIFT_DIR, matDVEDRAG_DIR, matDVESIDE_DIR);
% s_t = [1:valMAXTIME].*valDELTIME;
% valAR = (valSPAN.^2)./valAREA;
% CL_U2D = CL_U.*((valAR + 2)/valAR);
% plot(s_t, CL_U2D, '--bs', 'LineWidth', 1);

load('m5_dt0.2_wfix2.mat');
[~, CL_U, ~] = fcnDGAMMADT(1, valDELTIME, strATYPE, matINTCIRC, valDENSITY, valRPM, valDIAM, valAREA, valUINF, matLIFTFREE, matLIFTIND, matDRAGIND, matSIDEFREE, matSIDEIND, matDVELIFT_DIR, matDVEDRAG_DIR, matDVESIDE_DIR);
s_t = [1:valMAXTIME].*valDELTIME;
valAR = (valSPAN.^2)./valAREA;
CL_U2D = CL_U.*((valAR + 2)/valAR);
plot(s_t, CL_U2D, '-.rd', 'LineWidth', 1);
delt_time(2) = valDELTIME;
ratio(2) = valDELTIME/(1/5);

load('m5_dt0.05_wfix2.mat');
skip = 5;
[~, CL_U, ~] = fcnDGAMMADT(skip, valDELTIME, strATYPE, matINTCIRC, valDENSITY, valRPM, valDIAM, valAREA, valUINF, matLIFTFREE, matLIFTIND, matDRAGIND, matSIDEFREE, matSIDEIND, matDVELIFT_DIR, matDVEDRAG_DIR, matDVESIDE_DIR);
s_t = [1:skip:valMAXTIME].*valDELTIME;
valAR = (valSPAN.^2)./valAREA;
CL_U2D = CL_U.*((valAR + 2)/valAR);
plot(s_t, CL_U2D, '--m^', 'LineWidth', 1);
delt_time(3) = valDELTIME;
ratio(3) = (valDELTIME*skip)/(1/5);

hold off
xlabel('Time (s)');
ylabel('Lift Coefficient');

legend('Kussner Function', ['\Delta_T = ', num2str(delt_time(1)), 's, \Deltax_w/\Deltax_c = ', num2str(ratio(1))], ...
    ['\Delta_T = ', num2str(delt_time(2)), 's, \Deltax_w/\Deltax_c = ', num2str(ratio(2))], ...
    ['\Delta_T = ', num2str(delt_time(3)), 's, \Deltax_w/\Deltax_c = ', num2str(ratio(3))],'Location','SouthEast')


