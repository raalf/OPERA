function [CT, vecDVETHRUST] = fcnRFORCES(valTIMESTEP, matF_FS, matF_IF, matF_ID, valDENSITY, valDIAM, valRPM)

len = size(matF_FS,1);

vecDVETHRUST = dot(matF_FS(:,:,valTIMESTEP), repmat([0 0 1], len, 1), 2) ...
    + dot(matF_IF(:,:,valTIMESTEP), repmat([0 0 1], len, 1), 2) ...
    + dot(matF_ID(:,:,valTIMESTEP), repmat([0 0 1], len, 1), 2);

%% Output
CT = nansum(vecDVETHRUST)./(valDENSITY.*(pi.*((valDIAM/2).^2)).*(((valDIAM/2).*(valRPM.*(pi/30))).^2));

end