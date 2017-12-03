function [H00, H10, H01, H20, H11, H02] = fcnHINTEGRAL(a,h,l1,l2)
    
H00 = fcnH00(a,h,l2) - fcnH00(a,h,l1);
H10 = fcnH10(a,h,l2) - fcnH10(a,h,l1);
H01 = fcnH01(a,h,l2) - fcnH01(a,h,l1);
H20 = fcnH20(a,h,l2) - fcnH20(a,h,l1);
H11 = fcnH11(a,h,l2) - fcnH11(a,h,l1);
H02 = fcnH02(a,h,l2) - fcnH02(a,h,l1);
  
end

function out = fcnH00(a,h,l)

% out = (1./h).*atan((a.*l)./(a.^2 + h.^2 + h.*sqrt(l.^2 + a.^2 + h.^2)));
out = (1./h).*atan2((a.*l), (a.^2 + h.^2 + h.*sqrt(l.^2 + a.^2 + h.^2)));

end

function out = fcnH10(a,h,l)

out = (l./sqrt(l.^2 + a.^2)).*log(sqrt(l.^2 + a.^2 + h.^2) + sqrt(l.^2 + a.^2)) - log(l + sqrt(l.^2 + a.^2 + h.^2)) - ((l.*log(h))./(sqrt(l.^2 + a.^2)));

end

function out = fcnH01(a,h,l)

out = (a./sqrt(l.^2 + a.^2)).*(log(h) - log(sqrt(l.^2 + a.^2 + h.^2) + sqrt(l.^2 + a.^2)));

end

function out = fcnH20(a,h,l)

% out = ((a.*l.*(l.^2 + a.^2 + h.^2 - h.*sqrt(l.^2 + a.^2 + h.^2)))./((l.^2 + a.^2).*sqrt(l.^2 + a.^2 + h.^2))) - (h.*atan(l./a)) + h.*atan((h.*l)./(a.*sqrt(l.^2 + a.^2 + h.^2)));
out = ((a.*l.*(l.^2 + a.^2 + h.^2 - h.*sqrt(l.^2 + a.^2 + h.^2)))./((l.^2 + a.^2).*sqrt(l.^2 + a.^2 + h.^2))) - (h.*atan2(l,a)) + h.*atan2((h.*l),(a.*sqrt(l.^2 + a.^2 + h.^2)));

end

function out = fcnH11(a,h,l)

out = ((-a.^2).*(l.^2 + a.^2 + h.^2 - h.*sqrt(l.^2 + a.^2 + h.^2)))./((l.^2 + a.^2).*sqrt(l.^2 + a.^2 + h.^2));

end

function out = fcnH02(a,h,l)

% out = -h.*atan(l./a) + h.*atan((h.*l)./(a.*sqrt(l.^2 + a.^2 + h.^2)));
out = -h.*atan2(l, a) + h.*atan2((h.*l), (a.*sqrt(l.^2 + a.^2 + h.^2)));

end