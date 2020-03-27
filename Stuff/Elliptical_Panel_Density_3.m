clear

cd G:\GIT\opera

tmpstrSPACING = 'HALFCOSINE';
tmpvalALPHA = 4;
tmpAR = 7;

tmpvalMAXTIME = 0;
tmpvalDELTIME = 1;

xtcr = 1;
% N
inp = 2:3:17;
% inp = 2:3:8;
tmpvecM = 5;

for ii = 1:length(inp)

    elliptical_wing_o_matic(inp(ii), tmpvecM, tmpvalALPHA, tmpstrSPACING, 'NORMAL', tmpvalDELTIME, tmpvalMAXTIME, tmpAR, xtcr, []);
    OPERA_MAIN
    axis off
    fcnFIG2LATEX(hFig1, 'w1', [10 6])

    clearvars -except ax ii hFig2 tmpCL_N time tmpCDi_N tmpCL_M tmpCDi_M tmpE_N inp tmpvecN tmpvecM tmpstrSPACING tmpvalALPHA tmpAR tmpvalMAXTIME tmpvalDELTIME xtcr
end

% % % M 
% % inp = 1:2:12; 
% % % inp = 1:2:4;
% % tmpvecN = 11;
% % for ii = 1:length(inp)
% %     elliptical_wing_o_matic(tmpvecN, inp(ii), tmpvalALPHA, tmpstrSPACING, 'HALFCOSINE', tmpvalDELTIME, tmpvalMAXTIME, tmpAR, xtcr , []);
% %     OPERA_MAIN
% %     tmpCL_M(ii) = CL(end);
% %     tmpCDi_M(ii) = CDi(end);
% %     clearvars -except ii tmpCL_N tmpCDi_N tmpCL_M tmpCDi_M inp tmpvecN tmpvecM tmpstrSPACING tmpvalALPHA tmpAR tmpvalMAXTIME tmpvalDELTIME xtcr
% % end
% 
% save('Elliptical_Panel_Density_Deltime_3.mat')

%%
% load('Elliptical_Panel_Density_Deltime_3.mat');
% 
% nvals = inp.*5;
% 
% hFig1 = figure(1);
% clf(1);
% plot(inp.*2, tmpE_N,'-k', 'LineWidth',2);
% xlabel('Number of Spanwise Elements','FontSize',10)
% ylim([min(tmpE_N) 0.996])
% grid minor
% box off
% % axis tight
% ylabel('Span efficiency')
% 
% yyaxis right
% plot(inp.*2, time, '--b', 'LineWidth',2)
% ylabel('Runtime (s)')
% ax = gca;
% ax.YColor = 'k';
% legend('Span efficiency','Runtime','Location','SouthEast')
% fcnFIG2LATEX(hFig1, 'ellip', [5 5])

