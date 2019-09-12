function [CT, vecDVETHRUST] = fcnPFORCES(strATYPE, valTIMESTEP, vecDVELIFT, vecDVEDRAG, vecDVESIDE, matDVELIFT_DIR, matDVEDRAG_DIR, matDVESIDE_DIR, valDENSITY, valDIAM, valRPM)

vecDVETHRUST = dot(matDVELIFT_DIR.*vecDVELIFT, repmat([0 0 1], length(vecDVELIFT), 1), 2) ...
    + dot(matDVEDRAG_DIR.*vecDVEDRAG, repmat([0 0 1], length(vecDVEDRAG), 1), 2) ...
    + dot(matDVESIDE_DIR.*vecDVESIDE, repmat([0 0 1], length(vecDVESIDE), 1), 2);

%% Output
if strcmpi(strATYPE, 'PROPELLER')
    CT = nansum(vecDVETHRUST)./(valDENSITY.*((valRPM/60).^2).*(valDIAM.^4));
else
    CT = nansum(vecDVETHRUST)./(valDENSITY.*(pi.*((valDIAM/2).^2)).*(((valDIAM/2).*(valRPM.*(pi/30))).^2));
end
fprintf('Timestep: %d\t\tCT = %0.5f\n', valTIMESTEP, CT);

end