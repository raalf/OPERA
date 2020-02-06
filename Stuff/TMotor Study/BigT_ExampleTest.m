clear
clc

Alpha = 15;
RPM = 3000;
RotorName='T-Motor 18in';

MotorName = 'Scorpion';
ESCName = 'KDE';

Static_P = 29.76;
R = 0.2286;
A = pi*R^2;

Direction = 'CCW';
Trans = eye(6);
Dz = 6.95275;
Trans(4,2)= Dz;
Trans(5,1)= -Dz;

dataRate = 8000;

%%
encoderzero = getAIN(8,1500,10000);

%%
fcn_eWriteAddress(46180,0,0); %Set PID Off
fcn_eWriteAddress(46000,3,0); %Set Throttle to 0

%%
gotData = 0;
while gotData == 0
    [Tare] = getLabJackRawData( {'AIN0','AIN1','AIN2','AIN3','AIN4','AIN5'}, dataRate, 3, 0 );
    if sum(isnan(Tare)) == 0
        gotData = 1;
    end
end

Bias = mean(Tare);

%%
% fcn_eWriteAddress(46180,0,0); %Set PID Off

if exist('DATA','var') == 0
    DATA = struct;
end

load('CalData.mat');
Cal = FT20647;
c = clock;

lbf_N = 4.44822;
lbin_Nm = 0.112984;

for n=1:length(RPM)
    
    fcn_eWriteAddress(46180,0,0); %Set PID Off
    m = length(DATA) + 1;
    
    fprintf('Go to %i RPM\n ',RPM(n));
    
    
    
    
    % if the propeller is spinning too slow, start with THR = 20%
    % before turning PID to on to reduce overshoot
    %     if getRPM < 100
    %         fcn_eWriteAddress(46000,3,3.7);
    %         pause(1);
    %
    %     end
    
    fcn_eWriteAddress(46004,3,RPM(n)); % Setting target RPM
    fcn_eWriteAddress(46180,0,1); %Setting PID on
    
    disp('Hold for 20 seconds')
    buffsize = 1000;
    bufft = nan(1,buffsize);
    buffvel = nan(1,buffsize);
    buffrpm = nan(1,buffsize);
    tic
    pulsecount = 1;
    while toc<20 %plot the tunnel velocity while waiting
        bufft = [bufft(2:end) toc];
        vel = getTunnelSpeed( Static_P*3386.39, 48.2, 4.52 );
        rpm_temp = fcn_eReadAddress(46002,3)/pulsecount;
        buffvel = [buffvel(2:end) vel];
        buffrpm = [buffrpm(2:end) rpm_temp];
        
        figure(123)
        clf(123)
        
        subplot(2,1,1)
        plot(bufft,buffvel,'b','LineWidth',2);
        ylabel('Vel, m/s')
        axis tight
        grid on
        box on
        
        subplot(2,1,2)
        hold on
        plot(bufft,buffrpm,'b','LineWidth',2);
        plot(bufft,ones(size(bufft)).*RPM(n),'--r','LineWidth',2);
        axis tight
        ylabel('RPM')
        grid on
        box on
        
        drawnow
    end
    
    %             pause(10);
    
    Press = getAIN(6);
    Tempf = getAIN(7);
    Throttle = getThrottle;
    
    arsp_a = 48.2;
    arsp_b = 4.52;
    q_Pa = Press*arsp_a+arsp_b;
    Temp_f = Tempf*100;
    Temp_C = (Temp_f - 32) * 5/9;
    Temp_K = Temp_C + 273.15;
    rho = (Static_P* 3386.39)/287.058/Temp_K;
    V   = sqrt(2*q_Pa/rho);
    
    
    
    % read current RPM
    fcn_eReadAddress(46002,3)
    
    
    fcn_eWriteAddress(46180,0,0); %Setting PID off
    
    gotData = 0;
    while gotData == 0
        [RawData] = getLabJackRawData( {'AIN0','AIN1','AIN2','AIN3','AIN4','AIN5','AIN8'}, dataRate, 60*4, 0 );
        if sum(isnan(RawData)) == 0
            gotData = 1;
        end
    end
    
    fcn_eWriteAddress(46180,0,1); %Setting PID on
    
    
    % Write measurements and rawData into DATA structure
    DATA(m).Time = c;
    DATA(m).dataRate = dataRate;
    DATA(m).RotorName=RotorName;
    DATA(m).Direction=Direction;
    DATA(m).Alpha = Alpha;
    DATA(m).RPM = RPM(n);
    DATA(m).PressureDiff = Press;
    DATA(m).Tempf = Temp_f;
    DATA(m).TempC = Temp_C;
    DATA(m).TempK = Temp_K;
    DATA(m).rho = rho;
    DATA(m).Velocity = V;
    DATA(m).R   = R;
    
    DATA(m).encoderzero = encoderzero;
    
    DATA(m).Dz = Dz;
    DATA(m).Trans = Trans;
    DATA(m).Cal = Cal;
    DATA(m).Bias = Bias;
    DATA(m).RawData = RawData;
    
    
    % postprocessing rawData
    % Convert raw encoder data into rotor angle
    encoder = RawData(:,7);
    Angle = encoder/(max(encoder)-min(encoder))*360;
    unwrapped_Angle = rad2deg(unwrap(deg2rad(Angle)));
    unwrapped_Angle = unwrapped_Angle - DATA(m).encoderzero/(max(encoder)-min(encoder))*360;
    Angle = mod(unwrapped_Angle, 360);
    
    % Calculate actual mean rpm during test
    RPM_Actual = ((unwrapped_Angle(end) - unwrapped_Angle(1))./360)./((length(unwrapped_Angle).*(1/dataRate))./60);
    
    % Calculate omega using actual measured rpm
    Omega = (RPM_Actual*pi)/30;
    
    
    % convert raw F/T sensor voltages to forces/torques
    FT = (Trans*Cal*(RawData(:,1:6)'-Bias'))';
    
    %
    DATA(m).Throttle = Throttle;
    DATA(m).RPM_Actual = RPM_Actual;
    DATA(m).Omega = Omega;
    DATA(m).mu  = DATA(m).Velocity / (DATA(m).Omega * DATA(m).R);
    
    DATA(m).Angle = Angle;
    
    FxT = lbf_N*mean(FT(:,1));
    FyT = lbf_N*mean(FT(:,2));
    FzT = lbf_N*mean(FT(:,3));
    TxT = lbin_Nm*mean(FT(:,4));
    TyT = lbin_Nm*mean(FT(:,5));
    TzT = lbin_Nm*mean(FT(:,6));
    
    DATA(m).FxT_N = FxT;
    DATA(m).FyT_N = FyT;
    DATA(m).FzT_N = FzT;
    DATA(m).TxT_Nm = TzT;
    DATA(m).TyT_Nm = TyT;
    DATA(m).TzT_Nm = TzT;
    
    DATA(m).T_N = FzT;
    DATA(m).Q_Nm = TzT;
    DATA(m).Fx_N = -FyT;
    DATA(m).Mx_Nm = -TyT;
    DATA(m).Fy_N = FxT;
    DATA(m).My_Nm = TxT;
    
    Fcoeff = 1/(rho*A*((Omega*R)^2));
    Mcoeff = 1/(rho*A*((Omega^2)*(R^3)));
    CQ = -TzT*Mcoeff;
    
    DATA(m).CT = DATA(m).T_N*Fcoeff;
    DATA(m).CQ = DATA(m).Q_Nm*Mcoeff;
    DATA(m).CP = DATA(m).CQ;
    DATA(m).CFx = DATA(m).Fx_N*Fcoeff;
    DATA(m).CMx = DATA(m).Mx_Nm*Mcoeff;
    DATA(m).CFy = DATA(m).Fy_N*Fcoeff;
    DATA(m).CMy = DATA(m).My_Nm*Mcoeff;
    
end
close(123);
fcn_eWriteAddress(46004,3,RPM(1));
pause(3);

%%
setPIDOFF;
% setTHR(0);
fcn_eWriteAddress(46000,3,0); %Set Throttle to 0

%% Plot
% A(A(:,4) > 360,:) = [];

valDIAM = R*2;
valRPM = RPM;
valJ = mean(DATA(end).Velocity)/((valRPM.*(pi/30)).*(valDIAM/2));

vecPOS_TUNNEL_OG = DATA(end).Angle;
vecPOS_TUNNEL = unwrap(DATA(end).Angle);

vecDENSITY = DATA(end).rho;

% CT_tunnel = A(:,7)./(vecDENSITY.*(pi.*((valDIAM/2).^2)).*(((valDIAM/2).*(valRPM.*(pi/30))).^2));
% error_tunnel = (1/4)./(vecDENSITY.*(pi.*((valDIAM/2).^2)).*(((valDIAM/2).*(valRPM.*(pi/30))).^2));
CT_tunnel = lbf_N.*FT(:,3);
torque = lbin_Nm.*FT(:,6);

CT_tunnel = CT_tunnel./(rho.*(pi.*((valDIAM/2).^2)).*(((valDIAM/2).*(valRPM.*(pi/30))).^2));

hFig99 = figure(99);
clf(99);

subplot(3,1,1)
scatter(vecPOS_TUNNEL, CT_tunnel, 10, 'sk')
% yyaxis right
% scatter(vecPOS_TUNNEL, vecPOS_TUNNEL_OG, 10, 'm^')
% scatter(vecPOS_TUNNEL, RawData(:,7), 10, 'bo')
grid minor
box on
axis tight
ylabel('C_T')

subplot(3,1,2)
% [vecPOS_TUNNEL_OG, idx] = sort(vecPOS_TUNNEL_OG, 'ascend');
% scatter(vecPOS_TUNNEL_OG, CT_tunnel(idx), 10, 'sk')
scatter(vecPOS_TUNNEL_OG, CT_tunnel, 10, 'sk')

grid minor
box on
axis tight
ylabel('C_T');
xlabel('Azimuth Location (Degrees)')

subplot(3,1,3)
scatter(vecPOS_TUNNEL_OG, torque, 10, 'sk')
grid minor
box on
axis tight
ylabel('Torque (Nm)');
xlabel('Azimuth Location (Degrees)')

%%
TestDate = datestr(now,'yyyy-mm-dd');
files = dir('*.mat');
% matname = [TestDate,'/',MotorName, '_', RotorName, '_', num2str(Alpha), '_',num2str(DATA(end).Velocity),'.mat'];
matname = [TestDate,'/',datestr(datetime(DATA(end).Time)), '_', MotorName, '_', ESCName, '_', RotorName, '_RPM', num2str(RPM), '_Alpha', num2str(Alpha), '_', num2str(DATA(end).Velocity), '.mat']; 
matname = regexprep(matname, ':', '.');

count = 0;
for i = 1:size(files,1)
    if startsWith(files(i).name,matname)
        count = count + 1;
    end
end

if count > 0
    matname = [matname, '_', num2str(count + 1)];
end
save(matname)