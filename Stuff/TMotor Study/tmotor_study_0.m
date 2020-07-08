clc
clear

cd ./../../

% SEQ_JJ = [0.1 0.15 0.2 0.25];
% SEQ_JJ = [0.1 0.15 0.2 0.25];
SEQ_JJ = 0.05;
SEQ_ALPHA = SEQ_JJ.*0;

CT_sweep = nan(200, length(SEQ_JJ));
for jj = 1:length(SEQ_JJ)
%     figure(1);
%     set(gcf, 'Position', get(0, 'Screensize'));
    OPERA_MAIN
    close all
    filename = ['Stuff/TMotor Study/Alpha 0 Results/TMotor_Relaxed_J', num2str(SEQ_JJ(jj)), '.mat'];
    save(filename);
    CT_sweep(:,jj) = CT;
    clearvars -except SEQ_JJ SEQ_ALPHA jj CT_sweep
end

save('TMotor_sweep_0.mat')

clear

% SEQ_JJ = [0.1 0.15 0.2 0.25];
SEQ_JJ = 0.05;
SEQ_ALPHA = SEQ_JJ.*0 + 15;

CT_sweep = nan(200, length(SEQ_JJ));
for jj = 1:length(SEQ_JJ)
%     figure(1);
%     set(gcf, 'Position', get(0, 'Screensize'));
    OPERA_MAIN
    close all
    filename = ['Stuff/TMotor Study/Alpha 15 Results/TMotor_Relaxed_J', num2str(SEQ_JJ(jj)), '.mat'];
    save(filename);
    CT_sweep(:,jj) = CT;
    clearvars -except SEQ_JJ SEQ_ALPHA jj CT_sweep
end

save('TMotor_sweep_15.mat')