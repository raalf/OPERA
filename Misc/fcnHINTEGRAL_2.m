function [H00, H10, H01, H20, H11, H02] = fcnHINTEGRAL(a,h,l1,l2)
    
H00 = ((1./h).*atan2((a.*l2), (a.^2 + h.^2 + h.*sqrt(l2.^2 + a.^2 + h.^2)))) - ...
      ((1./h).*atan2((a.*l1), (a.^2 + h.^2 + h.*sqrt(l1.^2 + a.^2 + h.^2))));

H10 = ((l2./sqrt(l2.^2 + a.^2)).*log(sqrt(l2.^2 + a.^2 + h.^2) + sqrt(l2.^2 + a.^2)) - log(l2 + sqrt(l2.^2 + a.^2 + h.^2)) - ((l2.*log(h))./(sqrt(l2.^2 + a.^2)))) - ...
      ((l1./sqrt(l1.^2 + a.^2)).*log(sqrt(l1.^2 + a.^2 + h.^2) + sqrt(l1.^2 + a.^2)) - log(l1 + sqrt(l1.^2 + a.^2 + h.^2)) - ((l1.*log(h))./(sqrt(l1.^2 + a.^2))));

H01 = ((a./sqrt(l2.^2 + a.^2)).*(log(h) - log(sqrt(l2.^2 + a.^2 + h.^2) + sqrt(l2.^2 + a.^2)))) - ...
      ((a./sqrt(l1.^2 + a.^2)).*(log(h) - log(sqrt(l1.^2 + a.^2 + h.^2) + sqrt(l1.^2 + a.^2))));

H20 = (((a.*l2.*(l2.^2 + a.^2 + h.^2 - h.*sqrt(l2.^2 + a.^2 + h.^2)))./((l2.^2 + a.^2).*sqrt(l2.^2 + a.^2 + h.^2))) - h.*atan2(l2, a) + h.*atan2((h.*l2), (a.*sqrt(l2.^2 + a.^2 + h.^2)))) - ...
      (((a.*l1.*(l1.^2 + a.^2 + h.^2 - h.*sqrt(l1.^2 + a.^2 + h.^2)))./((l1.^2 + a.^2).*sqrt(l1.^2 + a.^2 + h.^2))) - h.*atan2(l1, a) + h.*atan2((h.*l1), (a.*sqrt(l1.^2 + a.^2 + h.^2))));

H11 = (((-a.^2).*(l2.^2 + a.^2 + h.^2 - h.*sqrt(l2.^2 + a.^2 + h.^2)))./((l2.^2 + a.^2).*sqrt(l2.^2 + a.^2 + h.^2))) - ...
      (((-a.^2).*(l1.^2 + a.^2 + h.^2 - h.*sqrt(l1.^2 + a.^2 + h.^2)))./((l1.^2 + a.^2).*sqrt(l1.^2 + a.^2 + h.^2)));
  
H02 = (-h.*atan2(l2, a) + h.*atan2((h.*l2), (a.*sqrt(l2.^2 + a.^2 + h.^2)))) - ...
      (-h.*atan2(l1, a) + h.*atan2((h.*l1), (a.*sqrt(l1.^2 + a.^2 + h.^2))));
  
end