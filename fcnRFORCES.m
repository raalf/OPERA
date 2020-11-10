function [CT, vecDVETHRUST, CMx, CMy, CMz] = fcnRFORCES(matF_FS, matF_IF, matF_ID, matF_NC, vecROTORAXIS, vecROTORHUB, matCENTER, valDENSITY, valDIAM, valRPM, matROTORTRANS)

len = size(matF_FS,1);

%% thrust
vecDVETHRUST = dot(matF_FS, repmat(vecROTORAXIS, len, 1), 2) ...
    + dot(matF_IF, repmat(vecROTORAXIS, len, 1), 2) ...
    + dot(matF_ID, repmat(vecROTORAXIS, len, 1), 2) ...
    + dot(matF_NC, repmat(vecROTORAXIS, len, 1), 2);

%% roll and pitch moments
r = matCENTER - vecROTORHUB;

vecDVEMOMENT = cross(r, matF_FS + matF_IF + matF_ID + matF_NC, 2);
CMx = nansum(dot(vecDVEMOMENT,  repmat(matROTORTRANS(:,:,1), len, 1), 2))./...
    (valDENSITY.*(pi.*((valDIAM/2).^2)).*((abs(valRPM).*(pi/30)).^2).*(valDIAM/2).^3);
CMy = nansum(dot(vecDVEMOMENT,  repmat(matROTORTRANS(:,:,2), len, 1), 2))./...
    (valDENSITY.*(pi.*((valDIAM/2).^2)).*((abs(valRPM).*(pi/30)).^2).*(valDIAM/2).^3);
CMz = nansum(dot(vecDVEMOMENT,  repmat(matROTORTRANS(:,:,3), len, 1), 2))./...
    (valDENSITY.*(pi.*((valDIAM/2).^2)).*((abs(valRPM).*(pi/30)).^2).*(valDIAM/2).^3);

%% Output
CT = nansum(vecDVETHRUST)./(valDENSITY.*(pi.*((valDIAM/2).^2)).*(((valDIAM/2).*(abs(valRPM).*(pi/30))).^2));

end