clc
clear

% cd ./../../

J = [0 0.03501408748 0.07002817496 0.1050422624 0.1400563499 0.1750704374];

CT_sweep = nan(250, length(J));
for jj = 1:length(J)
    OPERA_MAIN
    close all
    filename = ['Stuff/TMotor Study/TMotor_Relaxed_J', num2str(J(jj)), '.mat'];
    save(filename);
    CT_sweep(:,jj) = CT;
    clearvars -except J jj CT_sweep
end

