clc
clear

% cd ./../../

J = [0.0416 0.1113 0.2055 0.2753 0.3249];

CT_sweep = nan(200, length(J));
for jj = 1:length(J)
    figure(1);
    set(gcf, 'Position', get(0, 'Screensize'));
    OPERA_MAIN
    close all
    filename = ['Stuff/TMotor Study/Alpha 0 Results/TMotor_Relaxed_J', num2str(J(jj)), '.mat'];
    save(filename);
    CT_sweep(:,jj) = CT;
    clearvars -except J jj CT_sweep
end

