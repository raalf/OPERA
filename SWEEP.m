clc
clear

SEQ_ALPHA = [15 15 15 15 0 0 0 0];
SEQ_JJ = [0.0848 0.1109 0.1480 0.2113 0.1080 0.1441 0.2043 0.2889];

for jjj = 1:length(SEQ_ALPHA)
    
    OPERA_MAIN;
    matname = ['TMotor_Relaxed_J',num2str(SEQ_JJ(jjj)),'.mat'];
    save(matname);
    clearvars -except SEQ_ALPHA SEQ_JJ jjj
end