function [CT, vecDVETHRUST] = fcnRFORCES(matF_FS, matF_IF, matF_ID, vecROTORAXIS, valDENSITY, valDIAM, valRPM)

len = size(matF_FS,1);

vecDVETHRUST = dot(matF_FS, repmat(vecROTORAXIS, len, 1), 2) ...
    + dot(matF_IF, repmat(vecROTORAXIS, len, 1), 2) ...
    + dot(matF_ID, repmat(vecROTORAXIS, len, 1), 2);

%% Output
CT = nansum(vecDVETHRUST)./(valDENSITY.*(pi.*((valDIAM/2).^2)).*(((valDIAM/2).*(abs(valRPM).*(pi/30))).^2));

end