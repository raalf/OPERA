clear

% cd G:\GIT\opera
% 
% tmpstrSPACING = 'NORMAL';
% tmpvalALPHA = 4;
% tmpAR = 7;
% 
% tmpvalMAXTIME = 60;
% tmpvalDELTIME = 1;
% 
% % N
% inp = 2:15;
% tmpvecM = 3;
% for i = 1:length(inp)
%     elliptical_wing_o_matic(inp(i), tmpvecM, tmpvalALPHA, tmpstrSPACING, tmpvalDELTIME, tmpvalMAXTIME, tmpAR);
%     OPERA_MAIN
%     tmpCL_N(i) = CL(end);
%     clearvars -except i tmpCL_N inp tmpvecM tmpstrSPACING tmpvalALPHA tmpAR tmpvalMAXTIME tmpvalDELTIME
% end
% 
% % M 
% inp = 1:8;
% tmpvecN = 10;
% for i = 1:length(inp)
%     elliptical_wing_o_matic(tmpvecN, inp(i), tmpvalALPHA, tmpstrSPACING, tmpvalDELTIME, tmpvalMAXTIME, tmpAR);
%     OPERA_MAIN
%     tmpCL_M(i) = CL(end);
%     clearvars -except i tmpCL_N tmpCL_M inp tmpvecM tmpvecN tmpstrSPACING tmpvalALPHA tmpAR tmpvalMAXTIME tmpvalDELTIME
% end

load('Elliptical_Panel_Density.mat');
% load('FW_CL_v_N.mat')

hFig1 = figure(1);
clf(1);
plot((1:8).*2,tmpCL_M,'-^k');
hold on
plot((2:15).*2,tmpCL_N,'--sk');
hold off
legend('Increasing Chordwise Elements (20 Spanwise Elements)', 'Increasing Spanwise Elements (6 Chordwise Elements)','Location','SouthEast')
xlabel('Number of Variable Elements','FontSize',15);
ylabel('C_L','FontSize',15);
grid minor
set(gca,'YMinorTick','on')

