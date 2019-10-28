function infl_loc = fcnDBIND_FARFIELD(dvenum, P, a, nu_xi, nu_eta, epts, fpl)

%% D
% MXQ = 5
% PAGE 184-185

% a.)
x2 = epts(:,2,:,2);
x1 = epts(:,2,:,1);
A2 = x2;
A1 = x1;

D12 = A2 - A1;
D13 = A2.*x2 - A1.*x1;
D14 = (x2 + x1).*D13 - (x1.*x2.*D12);
D15 = (x2 + x1).*D14 - (x1.*x2.*D13);
D16 = (x2 + x1).*D15 - (x1.*x2.*D14);

% b.)
x2 = epts(:,1,:,2);
x1 = epts(:,1,:,1);
A2 = x2;
A1 = x1;

D21 = A2 - A1;
D31 = A2.*x2 - A1.*x1;
D41 = (x2 + x1).*D31 - (x1.*x2.*D21);
D51 = (x2 + x1).*D41 - (x1.*x2.*D31);
D61 = (x2 + x1).*D51 - (x1.*x2.*D41);

%% G

G11 = zeros(size(a));
G12 = G11; G13 = G11; G14 = G11; G15 = G11;
G21 = G11; G22 = G11; G23 = G11; G24 = G11;
G31 = G11; G32 = G11; G33 = G11; G41 = G11;
G42 = G11; G51 = G11;

% PAGE 137
% 6.)

% a.)
idx = abs(nu_eta) <= abs(nu_xi);

% (i)
N = 1;
G11(idx) = (1./(N.*nu_xi(idx))).*D12(idx);
N = 2;
G12(idx) = (1./(N.*nu_xi(idx))).*D13(idx);
N = 3;
G13(idx) = (1./(N.*nu_xi(idx))).*D14(idx);
N = 4;
G14(idx) = (1./(N.*nu_xi(idx))).*D15(idx);
N = 5;
G15(idx) = (1./(N.*nu_xi(idx))).*D16(idx);

% (ii)
G21(idx) = (-nu_eta(idx)./nu_xi(idx)).*G12(idx) + (a(idx)./nu_xi(idx)).*G11(idx);
G22(idx) = (-nu_eta(idx)./nu_xi(idx)).*G13(idx) + (a(idx)./nu_xi(idx)).*G12(idx);
G23(idx) = (-nu_eta(idx)./nu_xi(idx)).*G14(idx) + (a(idx)./nu_xi(idx)).*G13(idx);
G24(idx) = (-nu_eta(idx)./nu_xi(idx)).*G15(idx) + (a(idx)./nu_xi(idx)).*G14(idx);

G31(idx) = (-nu_eta(idx)./nu_xi(idx)).*G22(idx) + (a(idx)./nu_xi(idx)).*G21(idx);
G32(idx) = (-nu_eta(idx)./nu_xi(idx)).*G23(idx) + (a(idx)./nu_xi(idx)).*G22(idx);
G33(idx) = (-nu_eta(idx)./nu_xi(idx)).*G24(idx) + (a(idx)./nu_xi(idx)).*G23(idx);

G41(idx) = (-nu_eta(idx)./nu_xi(idx)).*G32(idx) + (a(idx)./nu_xi(idx)).*G31(idx);
G42(idx) = (-nu_eta(idx)./nu_xi(idx)).*G33(idx) + (a(idx)./nu_xi(idx)).*G32(idx);

G51(idx) = (-nu_eta(idx)./nu_xi(idx)).*G42(idx) + (a(idx)./nu_xi(idx)).*G41(idx);


% b.)
idx = abs(nu_xi) < abs(nu_eta);

% (i)
M = 1;
G11(idx) = (-1./(M.*nu_eta(idx))).*D21(idx);
M = 2;
G21(idx) = (-1./(M.*nu_eta(idx))).*D31(idx);
M = 3;
G31(idx) = (-1./(M.*nu_eta(idx))).*D41(idx);
M = 4;
G41(idx) = (-1./(M.*nu_eta(idx))).*D51(idx);
M = 5;
G51(idx) = (-1./(M.*nu_eta(idx))).*D61(idx);

% (ii)
G12(idx) = (-nu_xi(idx)./nu_eta(idx)).*G21(idx) + (a(idx)./nu_eta(idx)).*G11(idx);
G22(idx) = (-nu_xi(idx)./nu_eta(idx)).*G31(idx) + (a(idx)./nu_eta(idx)).*G21(idx);
G32(idx) = (-nu_xi(idx)./nu_eta(idx)).*G41(idx) + (a(idx)./nu_eta(idx)).*G31(idx);
G42(idx) = (-nu_xi(idx)./nu_eta(idx)).*G51(idx) + (a(idx)./nu_eta(idx)).*G41(idx);

G13(idx) = (-nu_xi(idx)./nu_eta(idx)).*G22(idx) + (a(idx)./nu_eta(idx)).*G12(idx);
G23(idx) = (-nu_xi(idx)./nu_eta(idx)).*G32(idx) + (a(idx)./nu_eta(idx)).*G22(idx);
G33(idx) = (-nu_xi(idx)./nu_eta(idx)).*G42(idx) + (a(idx)./nu_eta(idx)).*G32(idx);

G14(idx) = (-nu_xi(idx)./nu_eta(idx)).*G23(idx) + (a(idx)./nu_eta(idx)).*G13(idx);
G24(idx) = (-nu_xi(idx)./nu_eta(idx)).*G33(idx) + (a(idx)./nu_eta(idx)).*G23(idx);

G15(idx) = (-nu_xi(idx)./nu_eta(idx)).*G24(idx) + (a(idx)./nu_eta(idx)).*G14(idx);


%% C
% MXQ = 5

% PAGE 137

C11 = reshape((1./(1 + 1)).*sum(a.*G11, 3),1,1,[]);
C21 = reshape((1./(2 + 1)).*sum(a.*G21, 3),1,1,[]);
C31 = reshape((1./(3 + 1)).*sum(a.*G31, 3),1,1,[]);
C41 = reshape((1./(4 + 1)).*sum(a.*G41, 3),1,1,[]);
C51 = reshape((1./(5 + 1)).*sum(a.*G51, 3),1,1,[]);

C12 = reshape((1./(1 + 2)).*sum(a.*G12, 3),1,1,[]);
C22 = reshape((1./(2 + 2)).*sum(a.*G22, 3),1,1,[]);
C32 = reshape((1./(3 + 2)).*sum(a.*G32, 3),1,1,[]);
C42 = reshape((1./(4 + 2)).*sum(a.*G42, 3),1,1,[]);

C13 = reshape((1./(1 + 3)).*sum(a.*G13, 3),1,1,[]);
C23 = reshape((1./(2 + 3)).*sum(a.*G23, 3),1,1,[]);
C33 = reshape((1./(3 + 3)).*sum(a.*G33, 3),1,1,[]);

C14 = reshape((1./(1 + 4)).*sum(a.*G14, 3),1,1,[]);
C24 = reshape((1./(2 + 4)).*sum(a.*G24, 3),1,1,[]);

C15 = reshape((1./(1 + 5)).*sum(a.*G15, 3),1,1,[]);


%% J
len = length(C11);
P = reshape(P,1,1,[]);
Phat = reshape(fpl',3,1,[])./P;

% J11
E2 = [zeros(1,1,len); zeros(1,1,len); reshape(C11,1,1,[])];
E4 = [zeros(1,1,len) zeros(1,1,len) 0.5.*C21; ...
    zeros(1,1,len) zeros(1,1,len) 0.5.*C12; ...
    0.5.*C21 0.5.*C12 zeros(1,1,len)];
E4Phat = dot(E4, repmat(permute(Phat, [2 1 3]), 3, 1, 1), 2);
J11 = ((1./(P.^3)).*(E2 - 3.*dot(E2, Phat, 2).*Phat)) + ((1./(P.^4)).*(-15.*dot(Phat, E4Phat, 2).*Phat + 6.*E4Phat));

% J21
E2 = [zeros(1,1,len); zeros(1,1,len); reshape(C21,1,1,[])];
E4 = [zeros(1,1,len) zeros(1,1,len) 0.5.*C31; ...
    zeros(1,1,len) zeros(1,1,len) 0.5.*C22; ...
    0.5.*C31 0.5.*C22 zeros(1,1,len)];
E4Phat = dot(E4, repmat(permute(Phat, [2 1 3]), 3, 1, 1), 2);
J21 = ((1./(P.^3)).*(E2 - 3.*dot(E2, Phat, 2).*Phat)) + ((1./(P.^4)).*(-15.*dot(Phat, E4Phat, 2).*Phat + 6.*E4Phat));

% J22
E2 = [zeros(1,1,len); zeros(1,1,len); reshape(C22,1,1,[])];
E4 = [zeros(1,1,len) zeros(1,1,len) 0.5.*C32; ...
    zeros(1,1,len) zeros(1,1,len) 0.5.*C23; ...
    0.5.*C32 0.5.*C23 zeros(1,1,len)];
E4Phat = dot(E4, repmat(permute(Phat, [2 1 3]), 3, 1, 1), 2);
J22 = ((1./(P.^3)).*(E2 - 3.*dot(E2, Phat, 2).*Phat)) + ((1./(P.^4)).*(-15.*dot(Phat, E4Phat, 2).*Phat + 6.*E4Phat));

% J12
E2 = [zeros(1,1,len); zeros(1,1,len); reshape(C12,1,1,[])];
E4 = [zeros(1,1,len) zeros(1,1,len) 0.5.*C22; ...
    zeros(1,1,len) zeros(1,1,len) 0.5.*C13; ...
    0.5.*C22 0.5.*C13 zeros(1,1,len)];
E4Phat = dot(E4, repmat(permute(Phat, [2 1 3]), 3, 1, 1), 2);
J12 = ((1./(P.^3)).*(E2 - 3.*dot(E2, Phat, 2).*Phat)) + ((1./(P.^4)).*(-15.*dot(Phat, E4Phat, 2).*Phat + 6.*E4Phat));

% J31
E2 = [zeros(1,1,len); zeros(1,1,len); reshape(C31,1,1,[])];
E4 = [zeros(1,1,len) zeros(1,1,len) 0.5.*C41; ...
    zeros(1,1,len) zeros(1,1,len) 0.5.*C32; ...
    0.5.*C41 0.5.*C32 zeros(1,1,len)];
E4Phat = dot(E4, repmat(permute(Phat, [2 1 3]), 3, 1, 1), 2);
J31 = ((1./(P.^3)).*(E2 - 3.*dot(E2, Phat, 2).*Phat)) + ((1./(P.^4)).*(-15.*dot(Phat, E4Phat, 2).*Phat + 6.*E4Phat));

% J13
E2 = [zeros(1,1,len); zeros(1,1,len); reshape(C13,1,1,[])];
E4 = [zeros(1,1,len) zeros(1,1,len) 0.5.*C23; ...
    zeros(1,1,len) zeros(1,1,len) 0.5.*C14; ...
    0.5.*C23 0.5.*C14 zeros(1,1,len)];
E4Phat = dot(E4, repmat(permute(Phat, [2 1 3]), 3, 1, 1), 2);
J13 = ((1./(P.^3)).*(E2 - 3.*dot(E2, Phat, 2).*Phat)) + ((1./(P.^4)).*(-15.*dot(Phat, E4Phat, 2).*Phat + 6.*E4Phat));

%%

infl_loc = nan(3,6,length(dvenum));
infl_loc(:,1,:) = J13;
infl_loc(:,2,:) = J12;
infl_loc(:,3,:) = J31;
infl_loc(:,4,:) = J21;
infl_loc(:,5,:) = J22;
infl_loc(:,6,:) = J11;

end