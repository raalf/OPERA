clc
clear

cd ./../../

J = linspace(0.05, 0.3, 6)

CT_sweep = nan(160, length(J));
for jj = 1:length(J)
%     figure(1);
%     set(gcf, 'Position', get(0, 'Screensize'));
    OPERA_MAIN
    close all
    filename = ['Stuff/TMotor Study/Alpha 0 Results/New/TMotor_Relaxed_J', num2str(J(jj)), '.mat'];
    save(filename);
    CT_sweep(:,jj) = CT;
    clearvars -except J jj CT_sweep
end

