%
%***********************************************************************
%* LIFTING_LINE - A multi-lifting-line method for design and analysis  *
%*                of nonplanar wing-configurations                     *
%*                                                                     *
%* Version V2.5, January 2016                                          *
%*                                                                     *
%* Copyright (C) 1986 - 2016 Karl Heinz Horstmann                      *
%*                                                                     *
%***********************************************************************
%*                                                                     *
%* GENERAL INFORMATION                                                 *
%* -------------------                                                 *
%*                                                                     *
%* Sub-File: 'A23IND.f'                                                *
%*                                                                     *
%* This program is free software; you can redistribute it and/or       *
%* modify it under the terms of the GNU General Public License         *
%* as published by the Free Software Foundation; either version 2      *
%* of the License, or (at your option) any later version.              *
%*                                                                     *
%* This program is distributed in the hope that it will be useful,     *
%* but WITHOUT ANY WARRANTY; without even the implied warranty of      *
%* MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the       *
%* GNU General Public License for more details.                        *
%*                                                                     *
%* You should have received a copy of the GNU General Public License   *
%* along with this program; if not, write to the Free Software         *
%* Foundation, Inc., 59 Temple Place - Suite 330, Boston,              *
%* MA  02111-1307, USA.                                                *
%*                                                                     *
%* In case of questions or errors, please contact:                     *
%* Carsten Liersch                                                     *
%*  German Aerospace Center (DLR)                                      *
%*  Institute of Aerodynamics and Flow Technology                      *
%*  Transport Aircraft                                                 *
%*  Lilienthalplatz 7                                                  *
%*  d - 38108 Braunschweig / Germany                                   *
%*  Tel.  : +49 (0)531 / 295 - 2434                                    *
%*  Fax . : +49 (0)531 / 295 - 2320                                    *
%*  E-Mail: carsten.liersch@dlr.de                                     *
%*                                                                     *
%***********************************************************************

%***********************************************************************
%* LIFTING_LINE Sub-file 'A23IND.f'                                    *
function [bloc, cloc] = a23ind2(endpoints, phi, yaw, fpl, k)

EPS0 = 000001;
EPS_SMALL = 1.0d-25;
EPS_SMALL = 1e-2

hspan1 = abs(endpoints(:,1,2) - endpoints(:,1,1))./2;
hspan2 = abs(endpoints(:,2,2) - endpoints(:,2,1))./2;

elia = hspan1.*cos(yaw) + hspan2.*sin(yaw);

%-- PFEILWINKEL FI
tanfi = tan(phi);
atfi = abs(tanfi);

% Transformation from HDVE local to vortex sheet local
eta_translation = mean(endpoints,3);
fp_0 = fpl - eta_translation;
yia = fp_0(:,1).*cos(yaw) + fp_0(:,2).*-sin(yaw);
xia = fp_0(:,1).*sin(yaw) + fp_0(:,2).*cos(yaw);
zia = fp_0(:,3);
azia = abs(zia);

d = xia - yia.*tanfi;
ad = abs(d);

%--- V-WINKEL NUE
eps_ml = sqrt(((yia - elia) - (yia + elia)).^2)./1000;

len = length(zia);

%  FAKTOREN FUER W2
tj1 = yia - elia;
tj2 = yia + elia;
zamin = 03.*elia;

idx = (ad < eps_ml & azia < zamin & atfi > eps_ml);
%ML Verschiebung des betrachteten Punktes auf 3% hinter den Wirbel,
%ML damit nicht singulaer
xia(idx,1) = xia(idx) + 0.03.*elia(idx);
d(idx,1) = xia(idx) - yia(idx).*tanfi(idx);

a2 = 1 + tanfi.*tanfi;
b2 = d.*tanfi;
ab2 = abs(b2);
c2 = d.*d + zia.*zia;
n11 = tj1.*tj1 + zia.*zia;
n12 = tj2.*tj2 + zia.*zia;
n21 = sqrt(tj1.*tj1.*a2+2.0.*tj1.*b2+c2);
n22 = sqrt(tj2.*tj2.*a2+2.0.*tj2.*b2+c2);
epss = d.*d - zia.*zia.*tanfi.*tanfi;
ro = sqrt(epss.*epss+4.0.*zia.*zia.*tanfi.*tanfi.*d.*d);

awurz = ro./2.0 + epss./2.0;

idx = (awurz <= eps_ml.*eps_ml);
awurz(idx) = 0;

beta1 = -sqrt(awurz);
bwurz = ro./2.0 - epss./2.0;

idx = (bwurz <= eps_ml.*eps_ml);
bwurz(idx) = 0;

beta2 = -sqrt(bwurz);
dif2 = abs(zia.*b2) - eps_ml.*eps_ml;

idx = (dif2 > 0);
beta2(idx) = beta2(idx).*zia(idx).*b2(idx)./abs(zia(idx).*b2(idx));

idx = (ro > eps_ml.*eps_ml);
ga1(idx,1) =(a2(idx).*zia(idx).*beta2(idx) + b2(idx).*beta1(idx))./ro(idx);
ga2(idx,1) =(a2(idx).*zia(idx).*beta1(idx) - b2(idx).*beta2(idx))./ro(idx);
de1(idx,1) =(b2(idx).*zia(idx).*beta2(idx) + c2(idx).*beta1(idx))./ro(idx);
de2(idx,1) =(b2(idx).*zia(idx).*beta1(idx) - c2(idx).*beta2(idx))./ro(idx);

g211 = zeros(len,1);
g221 = zeros(len,1);
g212 = zeros(len,1);
g222 = zeros(len,1);
b21z1 = zeros(len,1) + 0.5;
b22z1 = -b21z1;
c21z1 = yia;
c22z1 = -c21z1;
b21y1 = zeros(len,1);
b22y1 = zeros(len,1);
c21y1 = zeros(len,1);
c22y1 = zeros(len,1);
b21z2 = zeros(len,1);
b22z2 = zeros(len,1);
c21z2 = zeros(len,1);
c22z2 = zeros(len,1);
b21y2 = ones(len,1);
b22y2 = ones(len,1);
c21y2 = yia.*2.0;
c22y2 = 2.0.*yia;
arg21 = a2.*tj2 + b2 + sqrt(a2).*n22;

% if(ad < eps_ml && azia < eps_ml)
idx = (ad < eps_ml & azia < eps_ml);

%     if(atfi >= eps_ml)
idx2 = (atfi >= eps_ml) & idx;
g211(idx2) = -log(((tj2(idx2).*tanfi(idx2) + n22(idx2)).^2.*tj1(idx2).*tj1(idx2))./((tj1(idx2).*tanfi(idx2) + n21(idx2)).^2.*tj2(idx2).*tj2(idx2)));
g221(idx2) = -log(((-a2(idx2)./tanfi(idx2).*tj1(idx2) - n21(idx2)).^2.*tj2(idx2).*tj2(idx2))./((-a2(idx2)./tanfi(idx2).*tj2(idx2)-n22(idx2)).^2.*tj1(idx2).*tj1(idx2)));

% elseif(ad >= eps_ml || atfi >= eps_ml)
idx = (ad >= eps_ml | atfi >= eps_ml);
arg11(idx,1) = ((ga1(idx).*tj2(idx) + de1(idx) - n22(idx)).^2 + (ga2(idx).*tj2(idx) + de2(idx)).^2).*(tj1(idx).*tj1(idx) + zia(idx).*zia(idx));
arg12(idx,1) = ((ga1(idx).*tj1(idx) + de1(idx) - n21(idx)).^2 + (ga2(idx).*tj1(idx) + de2(idx)).^2).*(tj2(idx).*tj2(idx) + zia(idx).*zia(idx));

%   if(arg12 > EPS_SMALL)
%   else
arg1(idx,1) = arg11(idx)./EPS_SMALL;
idx2 = (arg12 > EPS_SMALL) & idx;
arg1(idx2) = arg11(idx2)./arg12(idx2);

z1(idx,1) = ga2(idx).*tj2(idx) + de2(idx);
n1(idx,1) = ga1(idx).*tj2(idx) + de1(idx) - n22(idx);
z2(idx,1) = ga2(idx).*tj1(idx) + de2(idx);
n2(idx,1) = ga1(idx).*tj1(idx) + de1(idx) - n21(idx);

b21y1(idx,1) = -zia(idx).*tanfi(idx).*beta1(idx)./2.0./ro(idx);
b21y2(idx,1) = -zia(idx).*tanfi(idx).*beta2(idx)./ro(idx);
b22y1(idx,1) = -d(idx).*beta2(idx)./2.0./ro(idx);
b22y2(idx,1) = -d(idx).*beta1(idx)./ro(idx);
b21z1(idx,1) = -d(idx).*beta1(idx)./2.0./ro(idx);
b21z2(idx,1) = -d(idx).*beta2(idx)./ro(idx);
b22z1(idx,1) = zia(idx).*tanfi(idx).*beta2(idx)./2.0./ro(idx);
b22z2(idx,1) = zia(idx).*tanfi(idx).*beta1(idx)./ro(idx);
c21y1(idx,1) =(d(idx) - yia(idx).*tanfi(idx)).*zia(idx).*beta1(idx)./ro(idx);
c21y2(idx,1) =(d(idx) - yia(idx).*tanfi(idx)).*zia(idx).*beta2(idx).*2.0./ro(idx);

e(idx,1) = zia(idx).*zia(idx).*tanfi(idx) + yia(idx).*d(idx);

c22y1(idx,1) = -e(idx).*beta2(idx)./ro(idx);
c22y2(idx,1) = -2.0.*e(idx).*beta1(idx)./ro(idx);
c21z1(idx,1) = -e(idx).*beta1(idx)./ro(idx);
c21z2(idx,1) = -2.0.*e(idx).*beta2(idx)./ro(idx);
c22z1(idx,1) = -(d(idx) - yia(idx).*tanfi(idx)).*zia(idx).*beta2(idx)./ro(idx);
c22z2(idx,1) = -(d(idx) - yia(idx).*tanfi(idx)).*zia(idx).*beta1(idx).*2.0./ro(idx);

%     if(arg1 > EPS_SMALL)
%     else
g211(idx) = log(EPS_SMALL);
idx2 = (arg1 > EPS_SMALL) & idx;
g211(idx2) = log(arg1(idx2));

%     if(azia <= eps_ml)
g212(idx & azia <= eps_ml) = 0;

idx2 = (azia > eps_ml) & idx;

idx3 = (tj2.*tj2 < eps_ml.*eps_ml);
art1(idx2 & idx3) = pi./2.0.*zia(idx2 & idx3)./azia(idx2 & idx3);
idx4 = idx2 & ~idx3;
% art1(idx4) = atan2(zia(idx4), tj2(idx4));
art1(idx4) = atan(zia(idx4)./tj2(idx4));
art1(idx4 & zia > 0 & tj2 < 0) = art1(idx4 & zia > 0 & tj2 < 0) + pi;
art1(idx4 & zia < 0 & tj2 < 0) = art1(idx4 & zia < 0 & tj2 < 0) - pi;

idx3 = (tj1.*tj1 < eps_ml.*eps_ml);
art2(idx2 & idx3) = pi./2.0.*zia(idx2 & idx3)./azia(idx2 & idx3);
idx4 = idx2 & ~idx3;
% art2(idx4) = atan2(zia(idx4), tj1(idx4));
art2(idx4) = atan(zia(idx4)./tj1(idx4));
art2(idx4 & zia > 0 & tj1 < 0) = art2(idx4 & zia > 0 & tj1 < 0) + pi;
art2(idx4 & zia < 0 & tj1 < 0) = art2(idx4 & zia < 0 & tj1 < 0) - pi;

% art3(idx2) = atan2(z1(idx2), n1(idx2));
art3(idx2) = atan(z1(idx2)./n1(idx2));
art3(idx2 & n1 < 0) = art3(idx2 & n1 < 0) + pi;
art3(idx2 & z1 < 0 & n1 > 0) = art3(idx2 & z1 < 0 & n1 > 0) + 2.*pi;

% art4(idx2) = atan2(z2(idx2), n2(idx2));
art4(idx2) = atan(z2(idx2)./n2(idx2));
art4(idx2 & n2 < 0) = art4(idx2 & n2 < 0) + pi;
art4(idx2 & z2 < 0 & n2 > 0) = art4(idx2 & z2 < 0 & n2 > 0) + 2.*pi;

g212(idx2) = art1(idx2) - art2(idx2) + art3(idx2) - art4(idx2);

g221(idx) = log(1./EPS_SMALL);
g221(idx & arg1 > EPS_SMALL) = log(1./arg1(idx & arg1 > EPS_SMALL));

g222(idx) = g212(idx);

arg22 = ones(len,1);
idx = arg21 ~= 0;
arg22(idx) = a2(idx).*tj1(idx) + b2(idx) + sqrt(a2(idx)).*n21(idx);
idx = arg22 ~= 0;
arg2 = abs(arg21(idx)./arg22(idx));

g231 = n22 - n21;
g232 = log(arg2);
g24 = g232;

g25 = log(abs(EPS_SMALL./n11));

idx = (abs(n11) > EPS_SMALL) & (abs(n12) > EPS_SMALL);
g25(idx) = log(abs(n12(idx)./n11(idx)));

idx = abs(n11) <= EPS_SMALL;
g25(idx) = log(abs(n12(idx)./EPS_SMALL));

idx = (abs(n11) <= EPS_SMALL) & (abs(n12) <= EPS_SMALL);
g25(idx) = 0;

% g26 = atan2((tj2.*zia-tj1.*zia), (zia.*zia+tj1.*tj2));
g26 = atan((tj2.*zia-tj1.*zia)./(zia.*zia+tj1.*tj2));
x2y2 = tj1.*tj2./(zia.*zia);
x2 = tj2./zia;
idx = (x2y2 < -1 & x2 > 0) & abs(zia) >= EPS0;
g26(idx) = g26(idx) + pi;
idx = (x2y2 < -1 & x2 < 0) & abs(zia) >= EPS0;
g26(idx) = g26(idx) - pi;

g26(abs(zia) < EPS0) = 0;

g27 = tj2 - tj1;
b26y = -1.0;
b24z = -tanfi./sqrt(a2);
b25z = -0.5;
c24y = 2.0.*tanfi./sqrt(a2).*zia;
c25y = zia;
c26y = -2.0.*yia;
c23z1 = 2.0.*tanfi./a2;
c23z2 = -2.0.*tanfi.*b2./sqrt(a2.*a2.*a2);
c24z = 2.0.*(d-yia.*tanfi)./sqrt(a2);
c25z = -yia;
c26z = -2.0.*zia;
c27z = 2.0;

idx = ab2 < eps_ml & c2 < eps_ml.*eps_ml;

atj21(idx,1) = abs(tj2(idx)./tj1(idx));
g24(idx,1) = abs(log(atj21(idx)));

b2y = -(g211.*b21y1 + g212.*b21y2 + g221.*b22y1 + g222.*b22y2 + g26.*b26y);
c2y = -(g211.*c21y1 + g212.*c21y2 + g221.*c22y1 + g222.*c22y2 + g24.*c24y + g25.*c25y + g26.*c26y);
b2z = g211.*b21z1 + g212.*b21z2 + g221.*b22z1 + g222.*b22z2 +g24.*b24z + g25.*b25z;
c2z = g211.*c21z1 + g212.*c21z2 + g221.*c22z1 + g222.*c22z2 +g231.*c23z1 + g232.*c23z2 + g24.*c24z + g25.*c25z +g26.*c26z + g27.*c27z;

idx = azia < eps_ml;
b2y(idx) = 0;
c2y(idx) = 0;

idx = azia < eps_ml & ad < eps_ml;
b2z(idx) = 0;
c2z(idx) = 0;

% Begin Correction BZ2 and CZ2, 10/2014 Hor, implemented CML
axia = abs(xia);
idx = azia < eps_ml & atfi < eps_ml & axia < eps_ml;
arg(idx) = abs(tj2(idx)./tj1(idx));
b2z(idx) = -log(arg);
c2z(idx) = -(2.*yia(idx).*log(arg(idx))-2.*(tj2(idx)-tj1(idx)));
% End Correction BZ2 and CZ2, 10/2014 Hor, implemented CML

bloc = [zeros(length(b2z), 1) b2y b2z];
cloc = [zeros(length(c2z), 1) c2y c2z];

%%
tempb(:,2) = bloc(:,1).*cos(yaw) + bloc(:,2).*sin(yaw);
tempb(:,1) = bloc(:,1).*-sin(yaw) + bloc(:,2).*cos(yaw);
bloc(:,1:2) = tempb;

tempc(:,2) = cloc(:,1).*cos(yaw) + cloc(:,2).*sin(yaw);
tempc(:,1) = cloc(:,1).*-sin(yaw) + cloc(:,2).*cos(yaw);
cloc(:,1:2) = tempc;

end

