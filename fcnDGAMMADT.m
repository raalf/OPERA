function [CT_U, CL_U, matDGAMMADT] = fcnDGAMMADT(matINTCIRC, valDENSITY, valRPM, valDIAM, valAREA, valUINF, matLIFTFREE, matLIFTIND, matDRAGIND, matSIDEFREE, matSIDEIND, matDVELIFT_DIR, matDVEDRAG_DIR, matDVESIDE_DIR)

len = size(matLIFTFREE,1);
for i = 1:len
    if i <= len-2
        matDGAMMADT(i,:) = (-matINTCIRC(i+2,:) + 4.*matINTCIRC(i+1,:) - 3.*matINTCIRC(i,:))./(2.*valDELTIME);
    elseif i == len - 1
        matDGAMMADT(i,:) = ((matINTCIRC(i+1,:) - matINTCIRC(i,:))./valDELTIME);
    elseif i == len
        matDGAMMADT(i,:) = ((matINTCIRC(i,:) - matINTCIRC(i-1,:))./valDELTIME);
    end
end

tmpLIFTFREE = matLIFTFREE + matDGAMMADT;
tmpDVELIFT = tmpLIFTFREE + matLIFTIND;
tmpDVEDRAG = matDRAGIND;
tmpDVESIDE = matSIDEFREE + matSIDEIND;

for i = 1:size(matLIFTFREE,1)
    tmpDVETHRUST(i,:) = dot(matDVELIFT_DIR(:,:,i).*tmpDVELIFT(i,:)', repmat([0 0 1], length(tmpDVELIFT(i,:)), 1), 2) ...
        + dot(matDVEDRAG_DIR(:,:,i).*tmpDVEDRAG(i,:)', repmat([0 0 1], length(tmpDVEDRAG(i,:)), 1), 2) ...
        + dot(matDVESIDE_DIR(:,:,i).*tmpDVESIDE(i,:)', repmat([0 0 1], length(tmpDVESIDE(i,:)), 1), 2);
    
    if strcmpi(strATYPE, 'PROPELLER')
        CT_U(i,1) = nansum(tmpDVETHRUST(i,:))./(valDENSITY.*((valRPM/60).^2).*(valDIAM.^4));
    else
        CT_U(i,1) = nansum(tmpDVETHRUST(i,:))./(valDENSITY.*(pi.*((valDIAM/2).^2)).*(((valDIAM/2).*(valRPM.*(pi/30))).^2));
    end
    
end

CL_U = sum(tmpLIFTFREE,2)./(0.5.*valDENSITY.*valAREA.*(valUINF^2));


end



% len = size(matLIFTFREE,1);
% lambda = 0.5;
% matDGAMMADT(1,:) = matINTCIRC(1,:).*0;
% for i = 3:len
%     if i <= len-2
%         matDGAMMADT(i,:) = ((-matINTCIRC(i+2,:) + 8.*matINTCIRC(i+1,:) - 8.*matINTCIRC(i-1,:) + matINTCIRC(i-2,:))./(12.*valDELTIME));
%     elseif i == len - 1
%         matDGAMMADT(i,:) = ((matINTCIRC(i+1,:) - matINTCIRC(i-1,:))./(2.*valDELTIME));
%     elseif i == len
%         matDGAMMADT(i,:) = ((matINTCIRC(i,:) - matINTCIRC(i-1,:))./valDELTIME);
%     end
% end

% len = size(matLIFTFREE,1);
% lambda = 0.5;
% matDGAMMADT(1,:) = matINTCIRC(1,:).*0;
% for i = 2:len
%     matDGAMMADT(i,:) = ((matINTCIRC(i,:) - matINTCIRC(i-1,:))./valDELTIME).*lambda + (1 - lambda).*matDGAMMADT(i-1,:);
% end