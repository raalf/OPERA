clc
clear

load('FLAN_PAN.mat');

for ij = 1:size(ANGS,1)
    FLAP = ANGS(ij,2);
    PITCH = ANGS(ij,3);
    
    str = sprintf('Stuff/Columbia/Steady_Sweep_%d.mat', ij);
    
    OPERA_MAIN
    
    close all
    save(str)
    
    clearvars -except ANGS ij
end