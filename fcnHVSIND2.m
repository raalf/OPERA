function [aloc, bloc, cloc] = fcnHVSIND2(endpoints, limits, ELIA, phi, yaw, fpl, k)

%% This is a copy of Horstmann's vortex sheet induction function from a23ind.f
% It does not yet use singularity factor, and it has not been checked for accuracy.

% T.D.K. 2016-12-11 Rothwell Street, Aurora, Ontario, Canada. L4G-0V8

EPS_SMALL = 1e-14;
eps0 = 0.000001;

eta_translation = mean(endpoints,3);

fp_0 = fpl - eta_translation;

YIA = fp_0(:,1).*cos(yaw) + fp_0(:,2).*sin(yaw);
XIA = fp_0(:,1).*sin(yaw) + fp_0(:,2).*cos(yaw);
ZIA = fp_0(:,3);

AZIA = abs(ZIA);

% From Horstmann LIFTING_LINE, a23ind.f V-WINKEL NUE
% I believe this is the length of an unswept leading edge
% Horstmann uses zeta as well, but he does this before transformation
% to local ref. Here, zeta = 0 always. Length should be the same in
% either case.
EPS = sqrt(((YIA - ELIA) - (YIA + ELIA)).^2)/1000;

% From Horstmann LIFTING_LINE, a23ind.f FAKTOREN FUER W2
% I believe this is 3% of the half span (I think ELIA is
% the span of the vortex sheet, but it could be half span?)
ZAMIN = 0.03*(ELIA*2);

ZIAs = ZIA.*ZIA;

TANFI = tan(phi);
ATFI = abs(TANFI);

len = length(YIA(:,1));

G211 = zeros(len,1);
G222 = zeros(len,1);
G212 = zeros(len,1);
G221 = zeros(len,1);
G26 = zeros(len,1);
G24 = zeros(len,1);
G25 = zeros(len,1);
G231 = zeros(len,1);
G232 = zeros(len,1);
G27 = zeros(len,1);

D = (XIA - YIA.*TANFI);
AD = abs(D);

% From Horstmann LIFTING_LINE, a23ind.f FAKTOREN FUER W2 (150)
% If the point lies on the leading edge, we move it 3% back
idx200 = D < EPS & abs(ZIA) < ZAMIN & abs(TANFI) > EPS;
XIA(idx200) = XIA(idx200) + ZAMIN(idx200);
D(idx200) = XIA(idx200) - YIA(idx200).*TANFI(idx200);

% Eqn A2-12
A2 = 1 + (TANFI.^2);
B2 = D.*TANFI;
AB2 = abs(B2);
C2 = D.^2 + ZIAs;

TJ1 = YIA + (limits(:,1).*ELIA);
TJ2 = YIA + (limits(:,2).*ELIA);

TJ1s = TJ1.*TJ1;
TJ2s = TJ2.*TJ2;

N11 = TJ1s + ZIAs;
N12 = TJ2s + ZIAs;

N21 = sqrt((TJ1s).*A2 + 2.*TJ1.*B2 + C2);
N22 = sqrt((TJ2s).*A2 + 2.*TJ2.*B2 + C2);

% Eqn A2-5
EPSS = (D.^2) - (ZIAs).*(TANFI).^2;
RO = sqrt(EPSS.^2 + 4.*(ZIAs).*(B2.^2));

% Horstmann way
AWURZ = (RO./2) + (EPSS./2);
BETA1(AWURZ - (EPS.^2) <= 0, 1) = 0;
BETA1(AWURZ - (EPS.^2) > 0, 1) = -sqrt(AWURZ(AWURZ - (EPS.^2) > 0));

BWURZ = (RO./2) - (EPSS./2);
BETA2(BWURZ - (EPS.^2) <= 0, 1) = 0;
BETA2(BWURZ - (EPS.^2) > 0, 1) = -sqrt(BWURZ(BWURZ - (EPS.^2) > 0));

zetab2 = ZIA.*B2;
idx_B2 = (abs(zetab2) >= EPS.^2);
BETA2(idx_B2) = BETA2(idx_B2).*(zetab2(idx_B2))./abs(zetab2(idx_B2));

GA1 = (1./RO).*(A2.*BETA2.*ZIA + B2.*BETA1);
GA2 = (1./RO).*(A2.*BETA1.*ZIA - B2.*BETA2);
DE1 = (1./RO).*(B2.*BETA2.*ZIA + C2.*BETA1);
DE2 = (1./RO).*(B2.*BETA1.*ZIA - C2.*BETA2);

% Horstman a23ind.f 106 (ish)
idx205 = RO <= EPS.^2;
GA1(idx205) = 0;
GA2(idx205) = 0;
DE1(idx205) = 0;
DE2(idx205) = 0;

ARG11 = (GA1.*TJ2 + DE1 - N22).^2 + ((GA2.*TJ2 + DE2).^2).*(TJ1.*TJ1 + ZIAs);
ARG12 = (GA1.*TJ1 + DE1 - N21).^2 + ((GA2.*TJ1 + DE2).^2).*(TJ2.*TJ2 + ZIAs);

ARG1(ARG12 > EPS_SMALL, 1) = ARG11(ARG12 > EPS_SMALL)./ARG12(ARG12 > EPS_SMALL);
ARG1(ARG12 <= EPS_SMALL, 1) = ARG11(ARG12 <= EPS_SMALL)./EPS_SMALL;

% if any(find(~real(ARG12))) == 1
%     disp('Issue with h_arg12 in fcnHVSIND2');
% end

Z1 = GA2.*TJ2 + DE2;
N1 = GA1.*TJ2 + DE1 - N22;
Z2 = GA2.*TJ1 + DE2;
N2 = GA1.*TJ1 + DE1 - N21;

B21Y1 = -ZIA.*TANFI.*BETA1./2./RO;
B21Y2 = -ZIA.*TANFI.*BETA2./RO;
B22Y1 = -D.*BETA2./2./RO;
B22Y2 = -D.*BETA1./RO;
B21Z1 = -D.*BETA1./2./RO;
B21Z2 = -D.*BETA2./RO;
B22Z1 = ZIA.*TANFI.*BETA2./2./RO;
B22Z2 = ZIA.*TANFI.*BETA1./RO;
C21Y1 = (D - YIA.*TANFI).*ZIA.*BETA1./RO;
C21Y2 = (D - YIA.*TANFI).*ZIA.*BETA2.*2./RO;
E = ZIAs.*TANFI + YIA.*D;
C22Y1 = -E.*BETA2./RO;
C22Y2 = -2.*E.*BETA1./RO;
C21Z1 = -E.*BETA1./RO;
C21Z2 = -2.*E.*BETA2./RO;
C22Z1 = -(D - YIA.*TANFI).*ZIA.*BETA2./RO;
C22Z2 = -(D - YIA.*TANFI).*ZIA.*BETA1.*2./RO;

G211 = log(ARG1);
idx260 = ARG1 < EPS_SMALL;

if any(ARG1(idx260)) == 1
    disp('Issue with h_arg1 in fcnHVSIND2');
end

idx500 = AD < EPS & AZIA < EPS & ATFI > EPS;
G211(idx500) = -log(((TJ2(idx500).*TANFI(idx500) + N22(idx500)).^2.*TJ1(idx500).*TJ1(idx500))./((TJ1(idx500).*TANFI(idx500) + N21(idx500)).^2.*TJ2(idx500).*TJ2(idx500)));

idx501 = AD < EPS & ATFI < EPS;
G211(idx501) = zeros(length(nonzeros(idx501)),1);

ART1 = atan(ZIA./TJ2);
ART1(ZIA > 0 & TJ2 < 0) = ART1(ZIA > 0 & TJ2 < 0) + pi;
ART1(ZIA < 0 & TJ2 < 0) = ART1(ZIA < 0 & TJ2 < 0) - pi;
ART1(TJ2s < EPS.^2) = pi./2*ZIA(TJ2s < EPS.^2)./abs(ZIA(TJ2s < EPS.^2));

ART2 = atan(ZIA./TJ1);
ART2(ZIA > 0 & TJ1 < 0) = ART2(ZIA > 0 & TJ1 < 0) + pi;
ART2(ZIA < 0 & TJ1 < 0) = ART2(ZIA < 0 & TJ1 < 0) - pi;
ART2(TJ1s < EPS.^2) = pi./2*ZIA(TJ1s < EPS.^2)./abs(ZIA(TJ1s < EPS.^2));

ART3 = atan(Z1./N1);
ART3(N1 < 0) = ART3(N1 < 0) + pi;
ART3(Z1 < 0 & N1 > 0) = ART3(Z1 < 0 & N1 > 0) + 2.*pi;

ART4 = atan(Z2./N2);
ART4(N2 < 0) = ART4(N2 < 0) + pi;
ART4(Z2 < 0 & N2 > 0) = ART4(Z2 < 0 & N2 > 0) + 2.*pi;

G212 = ART1 - ART2 + ART3 - ART4;

G221 = log(1./ARG1);
G221(idx260) = repmat(log(1./EPS_SMALL),length(nonzeros(idx260)),1);
G221(idx500) = -log(((-A2(idx500)./TANFI(idx500).*TJ1(idx500) - N21(idx500)).^2.*TJ2(idx500).*TJ2(idx500))./((-A2(idx500)./TANFI(idx500).*TJ2(idx500) - N22(idx500)).^2.*TJ1(idx500).*TJ1(idx500)));
G221(idx501) = zeros(length(nonzeros(idx501)),1);

G222 = G212;

idx601 = AD < EPS & ATFI < EPS;
zero601 = zeros(length(nonzeros(idx601)),1);
len6 = length(zero601);

G212(idx601) = zero601;
G222(idx601) = zero601;
B21Z1(idx601) = zero601 + repmat(0.5, len6, 1);
B22Z1(idx601) = zero601 + repmat(-0.5, len6, 1);
C21Z1(idx601) = YIA(idx601);
C22Z1(idx601) = -YIA(idx601);
B21Y1(idx601) = zero601;
B22Y1(idx601) = zero601;
C21Y1(idx601) = zero601;
C22Y1(idx601) = zero601;
B21Z2(idx601) = zero601;
B22Z2(idx601) = zero601;
C21Z2(idx601) = zero601;
C22Z2(idx601) = zero601;
B21Y2(idx601) = ones(len6,1);
B22Y2(idx601) = ones(len6,1);
C21Y2(idx601) = 2.*YIA(idx601);
C22Y2(idx601) = 2.*YIA(idx601);

ARG21 = A2.*TJ2 + B2 + sqrt(A2).*N22;
ARG22 = A2.*TJ1 + B2 + sqrt(A2).*N21;

ARG22(ARG21 == 0) = 1;
ARG22(ARG22 == 0) = 1;

ARG2 = abs(ARG21./ARG22);

G231 = N22 - N21;
G232 = log(ARG2);
G24 = G232;

G25 = log(abs(N12./N11));
G25(N11 < EPS_SMALL) = log(abs(N12(N11 < EPS_SMALL)./EPS_SMALL));
G25(N12 < EPS_SMALL) = log(abs(EPS_SMALL./N11(N12 < EPS_SMALL)));
G25(N11 < EPS_SMALL & N12 < EPS_SMALL) = 0;

G26 = atan((TJ2.*ZIA - TJ1.*ZIA)./(ZIA.*ZIA + TJ1.*TJ2));
X2Y2 = TJ1.*TJ2./ZIAs;
X2 = TJ2./ZIA;

idx305 = X2Y2 < -1 & X2 > 0;
G26(idx305) = G26(idx305) + pi;

idx306 = X2Y2 < -1 & X2 < 0;
G26(idx306) = G26(idx306) - pi;

idx304 = abs(ZIA) < eps0;
G26(idx304) = 0;

G27 = TJ2 - TJ1;

B26Y = zeros(len,1) - 1;
B24Z = -TANFI./sqrt(A2);
B25Z = zeros(len,1) - 0.5;
C24Y = 2.*TANFI./sqrt(A2).*ZIA;
C25Y = ZIA;
C26Y = -2.*YIA;
C23Z1 = 2.*TANFI./A2;
C23Z2 = -2.*TANFI.*B2./sqrt(A2.*A2.*A2);
C24Z = 2.*(D - YIA.*TANFI)./sqrt(A2);
C25Z = -YIA;
C26Z = -2.*ZIA;
C27Z = zeros(len,1) + 2;

idx370 = AB2 < EPS & C2 < EPS.^2;
ATJ21(idx370) = abs(TJ2(idx370)./TJ1(idx370));
G24(idx370) = abs(log(ATJ21(idx370))); 

B2Y = -(G211.*B21Y1 + G212.*B21Y2 + G221.*B22Y1 + G222.*B22Y2 + G26.*B26Y);
C2Y = -(G211.*C21Y1 + G212.*C21Y2 + G221.*C22Y1 + G222.*C22Y2 + G24.*C24Y + G25.*C25Y + G26.*C26Y);

B2Z = G211.*B21Z1 + G212.*B21Z2 + G221.*B22Z1 + G222.*B22Z2 + G24.*B24Z + G25.*B25Z;
C2Z = G211.*C21Z1 + G212.*C21Z2 + G221.*C22Z1 + G222.*C22Z2 + G231.*C23Z1 + G232.*C23Z2 + G24.*C24Z + G25.*C25Z + G26.*C26Z + G27.*C27Z;

idx400 = AZIA < EPS & AD < EPS;
B2Z(idx400) = 0;
C2Z(idx400) = 0;

idx401 = AZIA < EPS;
B2Y(idx401) = 0;
C2Y(idx401) = 0;

idx402 = AZIA < EPS & ATFI < EPS & abs(XIA) < EPS;
ARG(idx402, 1) = abs(TJ2(idx402)./TJ1(idx402));
B2Z(idx402) = -log(ARG(idx402));
C2Z(idx402) = -(2.*YIA(idx402).*log(ARG(idx402)) - 2.*(TJ2(idx402) - TJ2(idx402)));

aloc = zeros(len,3);
bloc = [zeros(len,1) B2Y B2Z];
cloc = [zeros(len,1) C2Y C2Z];

%% Rotate 90 degrees to appropriate direction if needed

tempb(:,2) = bloc(:,1).*cos(yaw) + bloc(:,2).*sin(yaw);
tempb(:,1) = bloc(:,1).*sin(yaw) + bloc(:,2).*cos(yaw);
bloc(:,1:2) = tempb;

tempc(:,2) = cloc(:,1).*cos(yaw) + cloc(:,2).*sin(yaw);
tempc(:,1) = cloc(:,1).*sin(yaw) + cloc(:,2).*cos(yaw);
cloc(:,1:2) = tempc;

end

