function [matWAKEGEOM, WADJE, WELST, WVLST, WDVE, WNELE, WEATT, WEIDX, WELOC,...
    WPLEX, WDVECT, WALIGN, WVATT, WVNORM, WCENTER, matWCOEFF] = fcnCREATEWAKE(valTIMESTEP, matNEWWAKE, matWAKEGEOM, matCOEFF, valWSIZE, matWCOEFF, vecTE, matEATT)
% This function creates wake elements, and finds the fcnTRIANG values associated with the wake.
% All new values are calculated for the entire wake every timestep, as most of them will change
% as the wing moves and as the wake is relaxed. Currently unsure of the time penalty.

% INPUT:
%   valTIMESTEP - Current timestep number
%   matNEWWAKE - Points from the new wake, in an n x 3 x 3 matrix ready for triangulation
%   matWAKEGEOM - Past points from the wake, in an m x 3 x 3 matrix ready for triangulation
% OUTPUT:
%   matWAKEGEOM - Matrix of all wake points in an n+m x 3 x 3 matrix ready for triangulation at the next time step
%   All other values from fcnTRIANG

if valTIMESTEP <= 1
    matWAKEGEOM = matNEWWAKE;
    [~, WADJE, WELST, WVLST, WDVE, WNELE, WEATT, WEIDX, WELOC, WPLEX, WDVECT, WALIGN, WVATT, WVNORM, WCENTER] = fcnTRIANG('wake',matWAKEGEOM);
    matWCOEFF = repmat(matCOEFF(nonzeros(matEATT(vecTE,:)),:),2,1);
else
    matWAKEGEOM = cat(1, matWAKEGEOM, matNEWWAKE);
    [~, WADJE, WELST, WVLST, WDVE, WNELE, WEATT, WEIDX, WELOC, WPLEX, WDVECT, WALIGN, WVATT, WVNORM, WCENTER] = fcnTRIANG('wake',matWAKEGEOM);
    matWCOEFF = [matWCOEFF; repmat(matCOEFF(nonzeros(matEATT(vecTE,:)),:),2,1)];
end




