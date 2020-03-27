clc
clear

cd ./../../

J = linspace(0.05, 0.1405, 4);

CT_sweep = nan(160, length(J));
for jj = 1:length(J)
    OPERA_MAIN
    close all
    filename = ['Stuff/TMotor Study/Alpha 90 Results/TMotor_Relaxed_J', num2str(J(jj)), '.mat'];
    save(filename);
    CT_sweep(:,jj) = CT;
    clearvars -except J jj CT_sweep
end

