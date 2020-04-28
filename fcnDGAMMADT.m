function [CT_U, CL_U, matDGAMMADT, matF_NC] = fcnDGAMMADT(skip, valNELE, strATYPE, valDELTIME, matSPANDIR_ALL, matUINF_ALL, matF_FS, matF_IF, matF_ID, matINTCIRC, vecLIFT_DIR, vecSIDE_DIR, vecDRAG_DIR, valDENSITY, valAREA, valSPAN, valDIAM, valRPM)

matDGAMMADT = nan(size(matINTCIRC));
len = size(matINTCIRC,2);

%%
for i = 1:skip:len
    if i <= len - skip
        matDGAMMADT(:,i) = ((matINTCIRC(:,i+skip) - matINTCIRC(:,i))./(valDELTIME.*skip));
    else
        matDGAMMADT(:,i) = ((matINTCIRC(:,i) - matINTCIRC(:,i-skip))./(valDELTIME.*skip));
    end
end
matDGAMMADT(:,1) = matDGAMMADT(:,1).*0;

% for i = 1:skip:len
%     if i <= len-(2*skip)
%         matDGAMMADT(i,:) = (-matINTCIRC(i+(2*skip),:) + 4.*matINTCIRC(i+skip,:) - 3.*matINTCIRC(i,:))./(2.*valDELTIME.*skip);
%     elseif i > len-(2*skip) && i <= len - skip
%         matDGAMMADT(i,:) = ((matINTCIRC(i+skip,:) - matINTCIRC(i,:))./(valDELTIME.*skip));
%     else
%         matDGAMMADT(i,:) = ((matINTCIRC(i,:) - matINTCIRC(i-skip,:))./(valDELTIME.*skip));
%     end
% end

% for i = 2:skip:len
%     matDGAMMADT(i,:) = ((matINTCIRC(i,:) - matINTCIRC(i-skip,:))./(valDELTIME.*skip));
% end

%% Adding noncirculatory component to the freestream forces

matF_NC = permute(matDGAMMADT, [1 3 2]).*cross(matUINF_ALL./(sqrt(sum(matUINF_ALL.^2,2))), matSPANDIR_ALL, 2);

if strcmpi(strATYPE{1}, 'WING')
    CLf = nansum(dot(matF_FS + matF_NC, repmat(vecLIFT_DIR, valNELE, 1, len), 2))./(0.5.*valDENSITY.*valAREA);
    CLi = nansum(dot(matF_IF, repmat(vecLIFT_DIR, valNELE, 1, len), 2))./(0.5.*valDENSITY.*valAREA);
    CL_U = CLf(:) + CLi(:);
    CT_U = [];
else
    vecDVETHRUST = nansum(dot(matF_FS, repmat([0 0 1], valNELE, 1, len), 2) ...
        + dot(matF_NC, repmat([0 0 1], valNELE, 1, len), 2) ...
        + dot(matF_IF, repmat([0 0 1], valNELE, 1, len), 2) ...
        + dot(matF_ID, repmat([0 0 1], valNELE, 1, len), 2), 1);
    
    if strcmpi(strATYPE, 'PROPELLER')
        CT_U = vecDVETHRUST(:)./(valDENSITY.*((valRPM/60).^2).*(valDIAM.^4));
    else
        CT_U = vecDVETHRUST(:)./(valDENSITY.*(pi.*((valDIAM/2).^2)).*(((valDIAM/2).*(valRPM.*(pi/30))).^2));
    end
    CL_U = [];
end

end

