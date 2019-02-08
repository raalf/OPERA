clear
addpath('addaxis')

% cd G:\GIT\opera
% 
% tmpstrSPACING = 'COSINE';
% tmpvalALPHA = 4;
% tmpAR = 7;
% 
% tmpvalMAXTIME = 60;
% tmpvalDELTIME = 1;
% 
% xtcr = 1;
% % N
% inp = 2:3:17;
% tmpvecM = 8;
% for i = 1:length(inp)
%     elliptical_wing_o_matic(inp(i), tmpvecM, tmpvalALPHA, tmpstrSPACING, 'HALFCOSINE', tmpvalDELTIME, tmpvalMAXTIME, tmpAR, xtcr);
%     OPERA_MAIN
%     tmpCL_N(i) = CL(end);
%     tmpCDi_N(i) = CDi(end);
%     clearvars -except i tmpCL_N tmpCDi_N tmpCL_M tmpCDi_M inp tmpvecN tmpvecM tmpstrSPACING tmpvalALPHA tmpAR tmpvalMAXTIME tmpvalDELTIME xtcr
% end
% 
% % M 
% inp = 1:2:12;
% tmpvecN = 11;
% for i = 1:length(inp)
%     elliptical_wing_o_matic(tmpvecN, inp(i), tmpvalALPHA, tmpstrSPACING, 'HALFCOSINE', tmpvalDELTIME, tmpvalMAXTIME, tmpAR, xtcr );
%     OPERA_MAIN
%     tmpCL_M(i) = CL(end);
%     tmpCDi_M(i) = CDi(end);
%     clearvars -except i tmpCL_N tmpCDi_N tmpCL_M tmpCDi_M inp tmpvecN tmpvecM tmpstrSPACING tmpvalALPHA tmpAR tmpvalMAXTIME tmpvalDELTIME xtcr
% end
% 
% save('Elliptical_Panel_Density_Deltime_1.mat')

%%
load('Elliptical_Panel_Density_Deltime_1.mat');

nvals = (2:3:17);

hFig1 = figure(1);
clf(1);
plot(nvals.*2, tmpCDi_N,'-^k');
xlabel('Number of Spanwise Elements','FontSize',10)
hold on
addaxis(nvals.*2, (tmpCL_N.^2)./(pi.*tmpAR.*tmpCDi_N),'-.dk')
addaxis(nvals.*2, tmpCL_N,'--sk')
addaxislabel(1, 'C_D_i');
addaxislabel(2, 'e');
addaxislabel(3, 'C_L');
grid minor
box on
axis tight
legend('C_D_i','e','C_L','Location','East');
AX=findall(0,'type','axes'); 
set(AX(1),'fontsize',8)
set(AX(2),'fontsize',8)
set(AX(3),'fontsize',8)

dim = [.45 .3 .5 .3];
str = sprintf('alpha = 4\nx_t/c_r = 1\nAR = 7\nM = 8');
annotation('textbox',dim,'String',str,'FitBoxToText','on');

%%

nvals = (1:2:12);

hFig2 = figure(2);
clf(2);
plot(nvals, tmpCDi_M, '-^k');
xlabel('Number of Chordwise Elements','FontSize',10)
hold on
addaxis(nvals, (tmpCL_M.^2)./(pi.*tmpAR.*tmpCDi_M),'-.dk')
addaxis(nvals, tmpCL_M,'--sk')
addaxislabel(1, 'C_D_i');
addaxislabel(2, 'e');
addaxislabel(3, 'C_L');
grid minor
box on
axis tight
legend('C_D_i','e','C_L','Location','East');
AX=findall(0,'type','axes'); 
set(AX(1),'fontsize',8)
set(AX(2),'fontsize',8)
set(AX(3),'fontsize',8)

dim = [.45 .3 .5 .3];
str = sprintf('alpha = 4\nx_t/c_r = 1\nAR = 7\nN = 22');
annotation('textbox',dim,'String',str,'FitBoxToText','on');

