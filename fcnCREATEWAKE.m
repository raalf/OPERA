function [matWAKEGEOM, matWADJE, matWELST, matWVLST, matWDVE, valWNELE, matWEATT, matWEIDX, matWELOC,...
    matWPLEX, matWDVECT, matWVATT, matWVNORM, matWCENTER, matWCOEFF, matWROTANG] = fcnCREATEWAKE(valTIMESTEP, matNEWWAKE, matWAKEGEOM, matCOEFF, valWSIZE, ...
    vecTE, vecTEDVE, matPLEX, matELOC, matELST, matDVE, matWCOEFF)
% This function creates wake elements, and finds the fcnTRIANG values associated with the wake.
% All new values are calculated for the entire wake every timestep, as most of them will change
% as the wing moves and as the wake is relaxed. Currently unsure of the time penalty.

% INPUT:
%   valTIMESTEP - Current timestep number
%   matNEWWAKE - Points from the new wake, in an n x 3 x 3 matrix ready for triangulation
%   matWAKEGEOM - Past points from the wake, in an m x 3 x 3 matrix ready for triangulation
% OUTPUT:
%   matWAKEGEOM - Matrix of all wake points in an n+m x 3 x 3 matrix ready for triangulation at the next time step
%   All other values from fcnTRIANG

if valTIMESTEP <= 1
    matWAKEGEOM = matNEWWAKE;
    [~, matWADJE, matWELST, matWVLST, matWDVE, valWNELE, matWEATT, matWEIDX, matWELOC, matWPLEX, matWDVECT, matWVATT, matWVNORM, matWCENTER, matWROTANG] = fcnTRIANG(matWAKEGEOM);
    
    vecWLEDVE = [(valWNELE - 2*valWSIZE + 1):(valWNELE - valWSIZE)]'; % Post trailing edge row of wake HDVEs
    vecWLE = matWEIDX(vecWLEDVE,1);
    vecWTEDVE = [(valWNELE - valWSIZE + 1):valWNELE]';
    vecWTE = matWEIDX(vecWTEDVE,1);
    
    matNEWWAKECOEFF = fcnDWAKENEW(valWNELE, matPLEX, vecTEDVE, valWSIZE, matWPLEX, matELOC, vecTE, vecWLEDVE, matCOEFF, matWELOC, vecWLE, matDVE, matELST, matWDVE, matWELST, matWEATT, matWEIDX, vecWTEDVE, vecWTE);
    
    matWCOEFF = matNEWWAKECOEFF;
    
else
    matWAKEGEOM = cat(1, matWAKEGEOM, matNEWWAKE);
    [~, matWADJE, matWELST, matWVLST, matWDVE, valWNELE, matWEATT, matWEIDX, matWELOC, matWPLEX, matWDVECT, matWVATT, matWVNORM, matWCENTER, matWROTANG] = fcnTRIANG(matWAKEGEOM);
    
    vecWLEDVE = [(valWNELE - 2*valWSIZE + 1):(valWNELE - valWSIZE)]'; % Post trailing edge row of wake HDVEs
    vecWLE = matWEIDX(vecWLEDVE,1);
    vecWTEDVE = [(valWNELE - valWSIZE + 1):valWNELE]';
    vecWTE = matWEIDX(vecWTEDVE,1);
    
    matNEWWAKECOEFF = fcnDWAKENEW(valWNELE, matPLEX, vecTEDVE, valWSIZE, matWPLEX, matELOC, vecTE, vecWLEDVE, matCOEFF, matWELOC, vecWLE, matDVE, matELST, matWDVE, matWELST, matWEATT, matWEIDX, vecWTEDVE, vecWTE);
    matNEWWAKECOEFF = matNEWWAKECOEFF(end - valWSIZE*2 + 1:end,:);
%     matWCOEFF = repmat(matNEWWAKECOEFF,valTIMESTEP,1); % Steady
%     matWCOEFF(1:(end - valWSIZE*2),[3 4 5]) = matWCOEFF(1:(end - valWSIZE*2),[3 4 5]).*0
    matWCOEFF = [matWCOEFF; matNEWWAKECOEFF];
end




