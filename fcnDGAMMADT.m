function [CT_U, CL_U, matDGAMMADT, matDVENC] = fcnDGAMMADT(skip, valDELTIME, strATYPE, matINTCIRC, valDENSITY, valRPM, ...
    valDIAM, valAREA, valUINF, matLIFTFREE, matLIFTIND, matDRAGIND, matSIDEFREE, matSIDEIND, matDVELIFT_DIR, matDVEDRAG_DIR, matDVESIDE_DIR, matSPANDIR_ALL, matUINF_ALL, vecTE, vecTEDVE)

% matDGAMMADT = zeros(size(matLIFTFREE));
% len = size(matLIFTFREE,1);
% for i = 1:skip:len
%     if i <= len-(2*skip)
%         matDGAMMADT(i,:) = (-matINTCIRC(i+(2*skip),:) + 4.*matINTCIRC(i+skip,:) - 3.*matINTCIRC(i,:))./(2.*valDELTIME.*skip);
%     elseif i > len-(2*skip) && i <= len - skip
%         matDGAMMADT(i,:) = ((matINTCIRC(i+skip,:) - matINTCIRC(i,:))./(valDELTIME.*skip));
%     else
%         matDGAMMADT(i,:) = ((matINTCIRC(i,:) - matINTCIRC(i-skip,:))./(valDELTIME.*skip));
%     end
% end
% matDGAMMADT(1,:) = matDGAMMADT(1,:).*0;

if size(matINTCIRC,1) > 10
    for i = 1:size(matINTCIRC,2)
        
        matINTCIRC(:,i) = smooth(matINTCIRC(:,i),'rloess');
        
    end
end

matDGAMMADT = nan(size(matLIFTFREE));
len = size(matLIFTFREE,1);
for i = 1:skip:len
    if i <= len-(2*skip)
        matDGAMMADT(i,:) = (-matINTCIRC(i+(2*skip),:) + 4.*matINTCIRC(i+skip,:) - 3.*matINTCIRC(i,:))./(2.*valDELTIME.*skip);
    elseif i > len-(2*skip) && i <= len - skip
        matDGAMMADT(i,:) = ((matINTCIRC(i+skip,:) - matINTCIRC(i,:))./(valDELTIME.*skip));
    else
        matDGAMMADT(i,:) = ((matINTCIRC(i,:) - matINTCIRC(i-skip,:))./(valDELTIME.*skip));
    end
end
matDGAMMADT(1,:) = zeros(size(matDGAMMADT(1,:)));

% Adding noncirculatory component to the freestream forces
matDVENC = permute(matDGAMMADT, [2 3 1]).*cross(matUINF_ALL(vecTEDVE,:,:)./(sqrt(sum(matUINF_ALL(vecTEDVE,:,:).^2,2))), matSPANDIR_ALL(vecTEDVE,:,:), 2);
tmpDVELIFT = matLIFTFREE + matLIFTIND + permute(dot(matDVENC, matDVELIFT_DIR, 2), [3 1 2]);
tmpDVEDRAG = matDRAGIND;
tmpDVESIDE = matSIDEFREE + matSIDEIND + permute(dot(matDVENC, matDVESIDE_DIR, 2), [3 1 2]);

for i = 1:skip:size(matLIFTFREE,1)
    tmpDVETHRUST(i,:) = dot(matDVELIFT_DIR(:,:,i).*tmpDVELIFT(i,:)', repmat([0 0 1], length(tmpDVELIFT(i,:)), 1), 2) ...
        + dot(matDVEDRAG_DIR(:,:,i).*tmpDVEDRAG(i,:)', repmat([0 0 1], length(tmpDVEDRAG(i,:)), 1), 2) ...
        + dot(matDVESIDE_DIR(:,:,i).*tmpDVESIDE(i,:)', repmat([0 0 1], length(tmpDVESIDE(i,:)), 1), 2);
    
    if strcmpi(strATYPE, 'PROPELLER')
        CT_U(i,1) = nansum(tmpDVETHRUST(i,:))./(valDENSITY.*((valRPM/60).^2).*(valDIAM.^4));
    else
        CT_U(i,1) = nansum(tmpDVETHRUST(i,:))./(valDENSITY.*(pi.*((valDIAM/2).^2)).*(((valDIAM/2).*(valRPM.*(pi/30))).^2));
    end
    
end

CL_U = sum(tmpDVELIFT,2)./(0.5.*valDENSITY.*valAREA.*(valUINF^2));
CL_U = CL_U(~isnan(CL_U));

end

