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
plot(((s+2)/2), c_l, '-k');
box on
axis tight
grid minor

hold on
% disp('THIS IS WITH APPARENT MASS!!!!!!!!!!')
load('m5_dt1_wfix2.mat');
[~, CL_U, ~] = fcnDGAMMADT(valDELTIME, strATYPE, matINTCIRC, valDENSITY, valRPM, valDIAM, valAREA, valUINF, matLIFTFREE, matLIFTIND, matDRAGIND, matSIDEFREE, matSIDEIND, matDVELIFT_DIR, matDVEDRAG_DIR, matDVESIDE_DIR);
s_t = [1:valMAXTIME].*valDELTIME;
valAR = (valSPAN.^2)./valAREA;
CL_U2D = CL_U.*((valAR + 2)/valAR);
plot(s_t, CL_U2D, '-ok', 'LineWidth', 1);

load('m5_dt0.5_wfix2.mat');
[~, CL_U, ~] = fcnDGAMMADT(valDELTIME, strATYPE, matINTCIRC, valDENSITY, valRPM, valDIAM, valAREA, valUINF, matLIFTFREE, matLIFTIND, matDRAGIND, matSIDEFREE, matSIDEIND, matDVELIFT_DIR, matDVEDRAG_DIR, matDVESIDE_DIR);
s_t = [1:valMAXTIME].*valDELTIME;
valAR = (valSPAN.^2)./valAREA;
CL_U2D = CL_U.*((valAR + 2)/valAR);
plot(s_t, CL_U2D, '--bs', 'LineWidth', 1);

load('m5_dt0.25_wfix2.mat');
[~, CL_U, ~] = fcnDGAMMADT(valDELTIME, strATYPE, matINTCIRC, valDENSITY, valRPM, valDIAM, valAREA, valUINF, matLIFTFREE, matLIFTIND, matDRAGIND, matSIDEFREE, matSIDEIND, matDVELIFT_DIR, matDVEDRAG_DIR, matDVESIDE_DIR);
s_t = [1:valMAXTIME].*valDELTIME;
valAR = (valSPAN.^2)./valAREA;
CL_U2D = CL_U.*((valAR + 2)/valAR);
plot(s_t, CL_U2D, '-.rd', 'LineWidth', 1);

load('m5_dt0.2_wfix2.mat');
[~, CL_U, ~] = fcnDGAMMADT(valDELTIME, strATYPE, matINTCIRC, valDENSITY, valRPM, valDIAM, valAREA, valUINF, matLIFTFREE, matLIFTIND, matDRAGIND, matSIDEFREE, matSIDEIND, matDVELIFT_DIR, matDVEDRAG_DIR, matDVESIDE_DIR);
s_t = [1:valMAXTIME].*valDELTIME;
valAR = (valSPAN.^2)./valAREA;
CL_U2D = CL_U.*((valAR + 2)/valAR);
plot(s_t, CL_U2D, '--m^', 'LineWidth', 1);

hold off
xlabel('Time (s)');
ylabel('Lift Coefficient');

legend('Kussner Function', ['\Delta_T = 1s, \Deltax_w/\Deltax_c = ', num2str(1/(1/5))], ['\Delta_T = 0.5s, \Deltax_w/\Deltax_c = ', num2str(0.5/(1/5))], ['\Delta_T = 0.25s, \Deltax_w/\Deltax_c = ', num2str(0.25/(1/5))], ['\Delta_T = 0.05s, \Deltax_w/\Deltax_c = ', num2str(0.05/(1/5))],'Location','SouthEast')


