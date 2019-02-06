clear

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
load('Elliptical_Panel_Density_Deltime_1.mat')
hFig1 = figure(1);
clf(1);
yyaxis left
plot((1:2:12).*2,tmpCL_M,'-^k');
hold on
plot((2:3:17).*2,tmpCL_N,'--sk');
hold off
ylabel('C_L','FontSize',15);

yyaxis right
hold on
plot((1:2:12).*2,tmpCDi_M,'-.^k');
plot((2:3:17).*2,tmpCDi_N,':sk');
hold off

legend('Increasing Chordwise Elements (20 Spanwise Elements)', 'Increasing Spanwise Elements (5 Chordwise Elements)','Location','SouthEast')
xlabel('Number of Variable Elements','FontSize',15);
ylabel('C_D_i','FontSize',15);
grid minor
set(gca,'YMinorTick','on')

hFig2 = figure(2);
clf(2);
hold on
plot((1:2:12).*2,(tmpCL_M.^2)./(pi.*7.*tmpCDi_M),'-^k');
plot((2:3:17).*2,(tmpCL_N.^2)./(pi.*7.*tmpCDi_N),'--sk');
hold off
legend('Increasing Chordwise Elements (20 Spanwise Elements)', 'Increasing Spanwise Elements (5 Chordwise Elements)','Location','NorthEast')
xlabel('Number of Variable Elements','FontSize',15);
ylabel('Span Efficiency Factor','FontSize',15);
grid minor
set(gca,'YMinorTick','on')

