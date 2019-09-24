clc
clear

% cd ./../../

J = [0 0.0539 0.0783 0.1181 0.1404];

CT_sweep = nan(120, length(J));
for jj = 1:length(J)
    OPERA_MAIN
    close all
    filename = ['Stuff/TMotor Study/Alpha 90 Results/TMotor_Relaxed_J', num2str(J(jj)), '.mat'];
    save(filename);
    CT_sweep(:,jj) = CT;
    clearvars -except J jj CT_sweep
end

