function [n,bb,Ksiia,Etaia,Zetaia,Elia,Wksi1,Weta1,Wzet1,Weta2,Wzet2]=A23IND(n,bb,Ksiia,Etaia,Zetaia,Elia,Wksi1,Weta1,Wzet1,Weta2,Wzet2,varargin);



EPS0=0.000001;
EPS_SMALL=1.0d-25;


%------- Initialize PI-Parameter (only possible as variable):-----------
pi = 4.*atan(1.);

Wksi1 = 0.0;
Weta1 = 0.0;
Wzet1 = 0.0;
Weta2 = 0.0;
Wzet2 = 0.0;


ksij1 = h02r_1(j);
ksij2 = h02r_5(j);
etaj1 = h02r_2(j);
etaj2 = h02r_6(j);
zetaj1 = h02r_3(j);
zetaj2 = h02r_7(j);

%--- V-WINKEL NUE

detaj = etaj1 - etaj2;
dzetaj = zetaj1 - zetaj2;
wurz3 = sqrt(detaj.*detaj+dzetaj.*dzetaj);
snuej = dzetaj./wurz3;
cnuej = detaj./wurz3;
eps_ml = wurz3./1000.0;

%-- PFEILWINKEL FI

dksij = ksij1 - ksij2;
tanfi = dksij./wurz3;
atfi = abs(tanfi);

%--- ELEMENTURSPRUNG

ksi0 =(ksij1+ksij2)./2.0;
eta0 =(etaj1+etaj2)./2.0;
zeta0 =(zetaj1+zetaj2)./2.0;

%--- ORTSVEKTORTRANSFORMATION

yj1 =(etaj1-eta0).*cnuej +(zetaj1-zeta0).*snuej;
yj2 =(etaj2-eta0).*cnuej +(zetaj2-zeta0).*snuej;
xia = Ksiia - ksi0;
yia =(Etaia-eta0).*cnuej +(Zetaia-zeta0).*snuej;
zia = -(Etaia-eta0).*snuej +(Zetaia-zeta0).*cnuej;
azia = abs(zia);

%--- FAKTOREN FUER W1

a1 = 1 + tanfi.*tanfi;
b1 = -yia - xia.*tanfi;
c1 = xia.*xia + yia.*yia + zia.*zia;
d = xia - yia.*tanfi;
ad = abs(d);
r1 = sqrt(yj1.*yj1.*a1+2.0.*yj1.*b1+c1);
r2 = sqrt(yj2.*yj2.*a1+2.0.*yj2.*b1+c1);
if( azia>=eps_ml || ad>=eps_ml );
    adr1 = ad./r1;
    adr2 = ad./r2;
    
    if( adr1>=0.001 && adr2>=0.001 );
        goto 100;
    end;
    if( azia>eps_ml );
        goto 100;
    end;
end;
a1x = 0.0;
a1y = 0.0;
a1z = 0.0;
b1x = 0.0;
b1y = 0.0;
b1z = 0.0;
c1x = 0.0;
c1y = 0.0;
c1z = 0.0;
goto 150;
delta =(xia-yia.*tanfi).^2 + zia.*zia.*(1.+tanfi.*tanfi);
g11 =((a1.*yj2+b1)./(delta.*r2)-(a1.*yj1+b1)./(delta.*r1));
g12 =(-(b1.*yj2+c1)./(delta.*r2)+(b1.*yj1+c1)./(delta.*r1));
g131 =((2.0.*b1.*b1-a1.*c1).*yj2+b1.*c1)./(a1.*delta.*r2);
g132 =((2.0.*b1.*b1-a1.*c1).*yj1+b1.*c1)./(a1.*delta.*r1);
g133z = sqrt(a1).*r2 + a1.*yj2 + b1;
g133n = sqrt(a1).*r1 + a1.*yj1 + b1;
if((abs(g133z)<EPS0) ||(abs(g133n)<EPS0) );
    goto 50;
end;
g133 =(1.0./sqrt(a1.*a1.*a1)).*log(abs(g133z./g133n));
g13 = g131 - g132 + g133;

%
a1x = -g11.*zia;
a1y = g11.*zia.*tanfi;
a1z = g11.*d;
b1x = -g12.*zia;
b1y = g12.*zia.*tanfi;
b1z = g12.*d;
c1x = -g13.*zia;
c1y = g13.*zia.*tanfi;
c1z = g13.*d;
%
%  FAKTOREN FUER W2
%
tj1 = yia - yj1;
tj2 = yia - yj2;
zamin = 0.03.*Elia;
if( ad<eps_ml && azia<zamin && atfi>eps_ml );
    %ML Verschiebung des betrachteten Punktes auf 3% hinter den Wirbel,
    %ML damit nicht singulaer
    xia = xia + 0.03.*Elia;
    d = xia - yia.*tanfi;
end;
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
if( awurz<=eps_ml.*eps_ml );
    awurz = 0.0;
end;
beta1 = -sqrt(awurz);
bwurz = ro./2.0 - epss./2.0;
if( bwurz<=eps_ml.*eps_ml );
    bwurz = 0.0;
end;
beta2 = -sqrt(bwurz);
dif2 = abs(zia.*b2) - eps_ml.*eps_ml;
if( dif2>0.0 );
    beta2 = beta2.*zia.*b2./abs(zia.*b2);
end;
if( ro>eps_ml.*eps_ml );
    ga1 =(a2.*zia.*beta2+b2.*beta1)./ro;
    ga2 =(a2.*zia.*beta1-b2.*beta2)./ro;
    de1 =(b2.*zia.*beta2+c2.*beta1)./ro;
    de2 =(b2.*zia.*beta1-c2.*beta2)./ro;
end;
if( ad<eps_ml && azia<eps_ml );
    if( atfi>=eps_ml );
        g211 = -log(((tj2.*tanfi+n22).^2.*tj1.*tj1)./((tj1.*tanfi+n21).^2.*tj2.*tj2));
        g221 = -log(((-a2./tanfi.*tj1-n21).^2.*tj2.*tj2)./((-a2./tanfi.*tj2-n22).^2.*tj1.*tj1));
        goto 200;
    end;
elseif( ad>=eps_ml || atfi>=eps_ml ) ;
    
    arg11 =((ga1.*tj2+de1-n22).^2+(ga2.*tj2+de2).^2).*(tj1.*tj1+zia.*zia);
    arg12 =((ga1.*tj1+de1-n21).^2+(ga2.*tj1+de2).^2).*(tj2.*tj2+zia.*zia);
    
    if( arg12>EPS_SMALL );
        arg1 = arg11./arg12;
    else;
        arg1 = arg11./EPS_SMALL;
        
    end;
    
    z1 = ga2.*tj2 + de2;
    n1 = ga1.*tj2 + de1 - n22;
    z2 = ga2.*tj1 + de2;
    n2 = ga1.*tj1 + de1 - n21;
    b21y1 = -zia.*tanfi.*beta1./2.0./ro;
    b21y2 = -zia.*tanfi.*beta2./ro;
    b22y1 = -d.*beta2./2.0./ro;
    b22y2 = -d.*beta1./ro;
    b21z1 = -d.*beta1./2.0./ro;
    b21z2 = -d.*beta2./ro;
    b22z1 = zia.*tanfi.*beta2./2.0./ro;
    b22z2 = zia.*tanfi.*beta1./ro;
    c21y1 =(d-yia.*tanfi).*zia.*beta1./ro;
    c21y2 =(d-yia.*tanfi).*zia.*beta2.*2.0./ro;
    e = zia.*zia.*tanfi + yia.*d;
    c22y1 = -e.*beta2./ro;
    c22y2 = -2.0.*e.*beta1./ro;
    c21z1 = -e.*beta1./ro;
    c21z2 = -2.0.*e.*beta2./ro;
    c22z1 = -(d-yia.*tanfi).*zia.*beta2./ro;
    c22z2 = -(d-yia.*tanfi).*zia.*beta1.*2.0./ro;
    
    if( arg1>EPS_SMALL );
        g211 = log(arg1);
    else;
        g211 = log(EPS_SMALL);
        
    end;
    
    if( azia<=eps_ml );
        g212 = 0.0;
    else;
        if( tj2.*tj2<eps_ml.*eps_ml );
            art1 = pi./2.0.*zia./azia;
        else;
            art1 = atan(zia./tj2);
            if( zia>0.0 && tj2<0.0 );
                art1 = art1 + pi;
            end;
            if( zia<0.0 && tj2<0.0 );
                art1 = art1 - pi;
            end;
        end;
        if( tj1.*tj1<eps_ml.*eps_ml );
            art2 = pi./2.0.*azia./zia;
        else;
            art2 = atan(zia./tj1);
            if( zia>0.0 && tj1<0.0 );
                art2 = art2 + pi;
            end;
            if( zia<0.0 && tj1<0.0 );
                art2 = art2 - pi;
            end;
        end;
        art3 = atan(z1./n1);
        if( n1<0.0 );
            art3 = art3 + pi;
        end;
        if( z1<0.0 && n1>0.0 );
            art3 = art3 + 2.0.*pi;
        end;
        art4 = atan(z2./n2);
        if( n2<0.0 );
            art4 = art4 + pi;
        end;
        if( z2<0.0 && n2>0.0 );
            art4 = art4 + 2.0.*pi;
        end;
        g212 = art1 - art2 + art3 - art4;
    end;
    
    if( arg1>EPS_SMALL );
        g221 = log(1.0./arg1);
    else;
        g221 = log(1.0./EPS_SMALL);
        
    end;
    
    g222 = g212;
    goto 250;
end;
g211 = 0.0;
g221 = 0.0;
g212 = 0.0;
g222 = 0.0;
b21z1 = 0.5;
b22z1 = -0.5;
c21z1 = yia;
c22z1 = -yia;
b21y1 = 0.0;
b22y1 = 0.0;
c21y1 = 0.0;
c22y1 = 0.0;
b21z2 = 0.0;
b22z2 = 0.0;
c21z2 = 0.0;
c22z2 = 0.0;
b21y2 = 1.0;
b22y2 = 1.0;
c21y2 = yia.*2.0;
c22y2 = 2.0.*yia;
arg21 = a2.*tj2 + b2 + sqrt(a2).*n22;
if( arg21~=0 );
    arg22 = a2.*tj1 + b2 + sqrt(a2).*n21;
    if( arg22~=0 );
        arg2 = abs(arg21./arg22);
    else;
        arg2 = 1.0;
    end;
else;
    arg2 = 1.0;
end;
g231 = n22 - n21;
g232 = log(arg2);
g24 = g232;

if((abs(n11)>EPS_SMALL) &&(abs(n12)>EPS_SMALL) );
    g25 = log(abs(n12./n11));
else;
    
    if((abs(n11)<=EPS_SMALL) &&(abs(n12)<=EPS_SMALL) );
        g25 = 0.;
        
    elseif( abs(n11)<=EPS_SMALL ) ;
        g25 = log(abs(n12./EPS_SMALL));
        
    else;
        g25 = log(abs(EPS_SMALL./n11));
        
    end;
    
end;

if( abs(zia)<EPS0 );
    g26 = 0.0;
else;
    g26 = atan((tj2.*zia-tj1.*zia)./(zia.*zia+tj1.*tj2));
    x2y2 = tj1.*tj2./(zia.*zia);
    x2 = tj2./zia;
    if( x2y2<-1.0 && x2>0.0 );
        g26 = g26 + pi;
    end;
    if( x2y2<-1.0 && x2<0.0 );
        g26 = g26 - pi;
    end;
end;
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
if( ab2<eps_ml && c2<eps_ml.*eps_ml );
    atj21 = abs(tj2./tj1);
    g24 = log(atj21);
    g24 = abs(g24);
end;

b2y = -(g211.*b21y1+g212.*b21y2+g221.*b22y1+g222.*b22y2+g26.*b26y);
c2y = -(g211.*c21y1+g212.*c21y2+g221.*c22y1+g222.*c22y2+g24.*c24y+g25.*c25y+g26.*c26y);
b2z = g211.*b21z1 + g212.*b21z2 + g221.*b22z1 + g222.*b22z2 +g24.*b24z + g25.*b25z;
c2z = g211.*c21z1 + g212.*c21z2 + g221.*c22z1 + g222.*c22z2 +g231.*c23z1 + g232.*c23z2 + g24.*c24z + g25.*c25z +g26.*c26z + g27.*c27z;

if( azia<eps_ml && ad<eps_ml );
    a1z = 0.0;
    b1z = 0.0;
    c1z = 0.0;
    b2z = 0.0;
    c2z = 0.0;
end;
if( azia<eps_ml );
    a1y = 0.0;
    b1y = 0.0;
    c1y = 0.0;
    b2y = 0.0;
    c2y = 0.0;
end;

% Begin Correction BZ2 and CZ2, 10/2014 Hor, implemented CML
axia = abs(xia);
if( azia<eps_ml && atfi<eps_ml && axia<eps_ml );
    arg = abs(tj2./tj1);
    b2z = -log(arg);
    c2z = -(2.0.*yia.*log(arg)-2.0.*(tj2-tj1));
end;
% End Correction BZ2 and CZ2, 10/2014 Hor, implemented CML

aksi1 = a1x;
aeta1 = a1y.*cnuej - a1z.*snuej;
azet1 = a1y.*snuej + a1z.*cnuej;
bksi1 = b1x;
beta11 =(b1y).*cnuej -(b1z).*snuej;
bzet1 =(b1y).*snuej +(b1z).*cnuej;
cksi1 = c1x;
ceta1 =(c1y).*cnuej -(c1z).*snuej;
czet1 =(c1y).*snuej +(c1z).*cnuej;
beta22 = b2y.*cnuej - b2z.*snuej;
bzet2 = b2y.*snuej + b2z.*cnuej;
ceta2 = c2y.*cnuej - c2z.*snuej;
czet2 = c2y.*snuej + c2z.*cnuej;

azj = h04r_1(j);
bzj = h04r_2(j);
czj = h04r_3(j);
dwksi1 =(azj.*aksi1+bzj.*bksi1+czj.*cksi1);
dweta1 =(azj.*aeta1+bzj.*beta11+czj.*ceta1);
dwzet1 =(azj.*azet1+bzj.*bzet1+czj.*czet1);
dweta2 = bzj.*beta22 + czj.*ceta2;
dwzet2 = bzj.*bzet2 + czj.*czet2;
Wksi1 = Wksi1 + dwksi1;
Weta1 = Weta1 + dweta1;
Wzet1 = Wzet1 + dwzet1;
Weta2 = Weta2 + dweta2;
Wzet2 = Wzet2 + dwzet2;

Wksi1 = Wksi1.*bb./4.0./pi;
Weta1 = Weta1.*bb./4.0./pi;
Wzet1 = Wzet1.*bb./4.0./pi;
Weta2 = Weta2.*bb./4.0./pi;
Wzet2 = Wzet2.*bb./4.0./pi;



end