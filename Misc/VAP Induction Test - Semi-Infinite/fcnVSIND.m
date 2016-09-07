function [al, bl, cl] = fcnVSIND(endpoints, phi, eta_VS, fpl, k)

% INPUT:
%   endpoints - nx3x2 matrix, where (:,:,1) are first points and (:,:,2) are second. Maybe go around in clockwise direction
%   phi - nx1 vector of leading edge sweep, direction matters
%   eta_VS - nx3 vector of [1 0 0], eta of VS is same as eta of HDVE, or [0 1 0] where eta of VS is xsi of HDVE
%           In the local reference frame, vortex sheets will always go to the left, or down
%   fpl - nx3 matrix of field points in local coordinates of HDVE
%   k - singularity factor, not sure if we are using this now   
%
% OUTPUT:
%   al, bl, cl - local a, b, c for calculation of induced velocity

% Most equation labels are from Appendix 2 of Bramesfeld's thesis

% First we need to translate and rotate to the get the fp coordinates in the VS reference frame (_0)
eta_translation = mean(endpoints,3);

idx = ismember(eta_VS, [0 1 0], 'rows');

% half-span of the elements
hspan = endpoints(:,1,:);
hspan(idx,:,:) = endpoints(idx,2,:);
hspan = mean(hspan,3); % Do the negatives matter?

fp_0 = fpl - eta_translation;
fp_0(:,2) = -fp_0(:,2); % flipping xsi direction so it goes downwards in local ref frame 
fp_0(idx,:) = [-fp_0(idx,2) fp_0(idx,1) fp_0(idx,3)];

eta_0 = fp_0(:,1); 
xsi_0 = fp_0(:,2); 
zeta_0 = fp_0(:,3);


len = length(eta_0(:,1));

% Eqn A2-12
a2 = 1 + (tan(phi).^2);
b2 = (xsi_0 - eta_0.*tan(phi)).*tan(phi);
c2 = (xsi_0 - eta_0.*tan(phi)).^2 + zeta_0.^2;
t1 = eta_0 + hspan;
t2 = eta_0 - hspan;
rt_1 = sqrt((t1.^2).*a2 + 2.*t1.*b2 + c2);
rt_2 = sqrt((t2.^2).*a2 + 2.*t2.*b2 + c2);

% Eqn A2-5
eps = ((xsi_0 - eta_0.*tan(phi)).^2) - (zeta_0.^2).*(tan(phi)).^2;
rho = sqrt(eps.^2 + 4.*(zeta_0.^2).*(b2.^2));
beta1 = -sqrt((rho + eps)./2);
beta2 = -sqrt((rho - eps)./2);

% gamma1 = (1./rho).*(a2.*beta2.*zeta_0 + b2.*beta1);
% gamma2 = (1./rho).*(a2.*beta1.*zeta_0 - b2.*beta2);
% delta1 = (1./rho).*(b2.*beta2.*zeta_0 + c2.*beta1);
% delta2 = (1./rho).*(b2.*beta1.*zeta_0 - c2.*beta2);
% mu1_1 = ((gamma1.*t1 + delta1 - rt_1).^2 + (gamma2.*t1 + delta2).^2)./(k + t1.^2 + zeta_0.^2);
% mu1_2 = ((gamma1.*t2 + delta1 - rt_2).^2 + (gamma2.*t2 + delta2).^2)./(k + t2.^2 + zeta_0.^2);
% mu2_1 = atan2(zeta_0,t1) + atan2(gamma2.*t1 + delta2, gamma1.*t1 + delta1 - rt_1);
% mu2_2 = atan2(zeta_0,t2) + atan2(gamma2.*t2 + delta2, gamma1.*t2 + delta1 - rt_2);

% He changed the above equations to the below ones
gamma1 = (a2.*beta2.*zeta_0 + b2.*beta1);
gamma2 = (a2.*beta1.*zeta_0 - b2.*beta2);
delta1 = (b2.*beta2.*zeta_0 + c2.*beta1);
delta2 = (b2.*beta1.*zeta_0 - c2.*beta2);
mu1_1 = ((gamma1.*t1 + delta1 - rt_1.*rho).^2 + (gamma2.*t1 + delta2).^2);
mu1_2 = ((gamma1.*t2 + delta1 - rt_2.*rho).^2 + (gamma2.*t2 + delta2).^2);
mu2_1 = atan(zeta_0./t1) + atan((gamma2.*t1 + delta2)./(gamma1.*t1 + delta1 - rt_1.*rho));
mu2_2 = atan(zeta_0./t2) + atan((gamma2.*t2 + delta2)./(gamma1.*t2 + delta1 - rt_2.*rho));

% Eqn A2-3
G21 = ((beta1./(2.*rho)).*log(mu1_2) + (beta2./rho).*mu2_2) - ((beta1./(2.*rho)).*log(mu1_1) + (beta2./rho).*mu2_1);

% Eqn A2-4
G22 = ((1./xsi_0).*(-(beta2./(2.*rho)).*log(mu1_2) + (beta1./rho).*mu2_2)) ...
        - ((1./xsi_0).*(-(beta2./(2.*rho)).*log(mu1_1) + (beta1./rho).*mu2_1));
    
% Eqn A2-8
mu3_1 = a2.*t1 + b2 + sqrt(a2).*rt_1;
mu3_2 = a2.*t2 + b2 + sqrt(a2).*rt_2;

% Eqn A2-6
G23 = ((1./a2).*rt_2 - (b2./sqrt(a2.^3)).*log(mu3_2)) - ((1./a2).*rt_1 - (b2./sqrt(a2.^3)).*log(mu3_1));

% Eqn A2-7
G24 = ((1./sqrt(a2)).*log(mu3_2)) - ((1./sqrt(a2)).*log(mu3_1));

% Eqn A2-9
% G25 = (0.5.*log(k + t2.^2 + zeta_0.^2)) - (0.5.*log(k + t1.^2 + zeta_0.^2));
G25 = (0.5.*log(t2.^2 + zeta_0.^2)) - (0.5.*log(t1.^2 + zeta_0.^2));
% Eqn A2-10
G26 = ((1./zeta_0).*atan2(t2,zeta_0)) - (1./zeta_0).*atan2(t1,zeta_0);

% Eqn A2-11
G27 = t2 - t1;

% Eqn A2-13
b21 = -(xsi_0 - (eta_0.*tan(phi)));
b22 = (zeta_0.^2).*tan(phi);
b23 = zeros(len,1);
b24 = -tan(phi);
b25 = -ones(len,1);
b26 = zeros(len,1);
b27 = zeros(len,1);

c21 = -2.*((zeta_0.^2).*tan(phi) + eta_0.*(xsi_0 - eta_0.*tan(phi)));
c22 = -2.*(zeta_0.^2).*(xsi_0 - 2.*eta_0.*tan(phi)); % NOPE
c23 = 2.*tan(phi);
c24 = 2.*(xsi_0 - 2.*eta_0.*tan(phi));
c25 = -2.*eta_0;
c26 = -2.*(zeta_0.^2);
c27 = repmat(2,len,1);

% Eqn A2-2
b2_eta = -zeta_0.*(G21.*b24 + G22.*b21 + G26.*b25);
b2_zeta = G21.*b21 + G22.*b22 + G23.*b23 + G24.*b24 + G25.*b25 + G26.*b26 + G27.*b27;

c2_eta = -zeta_0.*(G21.*c24 + G22.*c21 + G24.*c23 + G25.*c27 + G26.*c25);
c2_zeta = G21.*c21 + G22.*c22 + G23.*c23 + G24.*c24 + G25.*c25 + G26.*c26 + G27.*c27;
% G 1 2 5

b2_xsi = zeros(size(b2_eta));
c2_xsi = zeros(size(c2_eta));

% a, b, c in local ref frame
bl = [b2_eta b2_xsi b2_zeta];
cl = [c2_eta c2_xsi c2_zeta];
al = zeros(size(bl));

end

