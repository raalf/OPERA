clc
clear


%% Roll
% 0.5, 1, 5 rad/s?

rot_type = 'ROLL';
SEQ_RATES = [-1 -0.5 -0.2 0.2 0.5 1];
% SEQ_VS = [10 15];
SEQ_VS = 10;
% SEQ_ALPHAS = [0 15];
SEQ_ALPHAS = 0;


for jjj = 1:length(SEQ_ALPHAS)
    for jjjj = 1:length(SEQ_VS)
        for jjjjj = 1:length(SEQ_RATES)
            
            OPERA_SUPREME_MAIN;
            matname = ['Runs/', rot_type, '_A', num2str(SEQ_ALPHAS(jjj)), '_V', num2str(SEQ_VS(jjjj)), '_R', num2str(SEQ_RATES(jjjjj)), '.mat'];
            save(matname);
            clearvars -except jjj jjjj jjjjj SEQ_ALPHAS SEQ_VS SEQ_RATES rot_type
            
        end
    end
end