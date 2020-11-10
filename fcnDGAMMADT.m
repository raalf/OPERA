function [matF_NC] = fcnDGAMMADT(skip, valTIMESTEP, valNELE, valDENSITY, valDELTIME, matINTCIRC, matDVECT)

matDGAMMADT = zeros(size(matINTCIRC));
len = size(matINTCIRC,2);

%%
for i = 2:skip:len
    if i <= len - skip
        matDGAMMADT(:,i) = ((matINTCIRC(:,i+skip) - matINTCIRC(:,i))./(valDELTIME.*skip));
    else
        matDGAMMADT(:,i) = ((matINTCIRC(:,i) - matINTCIRC(:,i-skip))./(valDELTIME.*skip));
    end
end

if valTIMESTEP > 1*skip
    matF_NC = valDENSITY.*permute(matDGAMMADT, [1 3 2]).*matDVECT(:,:,3);
else
    matF_NC = zeros(valNELE,3);
end

end

