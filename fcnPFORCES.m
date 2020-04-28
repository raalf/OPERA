function [CT, vecDVETHRUST] = fcnPFORCES(flagPRINT, strATYPE, valTIMESTEP, matF_FS, matF_IF, matF_ID, valDENSITY, valDIAM, valRPM, vecTSITER, vecRSQUARED)

% vecDVETHRUST = dot(matDVELIFT_DIR(:,:,valTIMESTEP).*vecDVELIFT', repmat([0 0 1], length(vecDVELIFT), 1), 2) ...
%     + dot(matDVEDRAG_DIR(:,:,valTIMESTEP).*vecDVEDRAG', repmat([0 0 1], length(vecDVEDRAG), 1), 2) ...
%     + dot(matDVESIDE_DIR(:,:,valTIMESTEP).*vecDVESIDE', repmat([0 0 1], length(vecDVESIDE), 1), 2);

len = size(matF_FS,1);

vecDVETHRUST = dot(matF_FS(:,:,valTIMESTEP), repmat([0 0 1], len, 1), 2) ...
    + dot(matF_IF(:,:,valTIMESTEP), repmat([0 0 1], len, 1), 2) ...
    + dot(matF_ID(:,:,valTIMESTEP), repmat([0 0 1], len, 1), 2);

%% Output
if strcmpi(strATYPE, 'PROPELLER')
    CT = nansum(vecDVETHRUST)./(valDENSITY.*((valRPM/60).^2).*(valDIAM.^4));
else
    CT = nansum(vecDVETHRUST)./(valDENSITY.*(pi.*((valDIAM/2).^2)).*(((valDIAM/2).*(valRPM.*(pi/30))).^2));
end

if flagPRINT
    fprintf('Timestep: %d\t\tCT = %0.5f\t\tIterations = %d\t\tR^2 = (%g, %g)\n', valTIMESTEP, CT, vecTSITER(valTIMESTEP), vecRSQUARED(valTIMESTEP,1), vecRSQUARED(valTIMESTEP,2));
end
end