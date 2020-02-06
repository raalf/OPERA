
function [Pxy,COH, F, Navg]=CalculatePSDandCohs(InputData, FFTParam)
%
%function [Pxy,  F]=calcpsdnicolet(InputData, FFTParam)
% This function calculates the psd of the data contained in InputData
% and determines the data format.  It recognizes one format:
% 1.  The MATLAB file exported by the Nicolet Odyssey DAQ System
%
% use :  [Pxy, Cxy, F]=calcpsdnicolet(InputData, FFTParam)
%
%	Inputs:
%		InputData  - name of data stored in the workspace STRING
%		fftparam - parameter structure for the ffts, contains   STRUCTURE
%				FFTLength - Length of the FFT
%				SampRate - Data acquisition sample rate
%				WinType - FFT Window Type (ie hanning)
%				PercOverlap - percent overlap for the FFT averaging
%				DFLAG - Data detrending flag (set to 'none')
%
%	Outputs:
%		Pxy - Auto & Cross Spectral Values for each channel pair - MATRIX
%		F - Frequency Values - VECTOR
%

% check if file name has been specified at all
% and if not set = to *

error = 0;

%
% determine the number of columns in InputData
%
NumChan = size(InputData,2);
NumSamp = size(InputData,1);
%
% Calculate the psd using MATLAB psd routine
%
if nargin==1
    
    % pop a GUI asking for the FFT parameters if not specified
    prompt={'Enter the FFT Length:',...
        'Enter the Sample Rate:',...
        'Enter the Window Type:',...
        'Enter the Number of Overlap Points:', ...
        'Enter the Data Detrend Mode:'};
    def={'8192','500000','hanning','0','none'};
    title='Input for Power Spectral Density Calculations';
    lineNo=1;
    answer=inputdlg(prompt,title,lineNo,def);
    FFTLength = str2num(char(answer{1}));
    SampRate = str2num(char(answer{2}));
    WinType = char(answer{3});
    PercOverlap = str2num(char(answer{4}));
    DFlag = char(answer{4});
else
    FFTLength = FFTParam.FFTLength;
    SampRate = FFTParam.SampRate;
    WinType = char(FFTParam.WinType);
    PercOverlap = FFTParam.PercOverlap;
    DFlag = char(FFTParam.DFlag);
end
Navg=floor((100/PercOverlap)*size(InputData,1)/FFTLength-1);
%
%Cxy = zeros(FFTLength/2+1, NumChan, NumChan);
Pxy = zeros(FFTLength/2+1, NumChan, NumChan);
COH= zeros(FFTLength/2+1, NumChan, NumChan);
%
window = eval([WinType,'(',num2str(FFTLength),')']);
WinCorr = norm(window)^2/sum(window)^2;
%
% loop through the channels and calculate PSD, CSD, and Coherence Functions
%
NOverlap=FFTLength*PercOverlap/100;
for ichan = 1:NumChan
    disp(['Calculating PSD for Channel ',num2str(ichan),' of ',num2str(NumChan),': Please be patient.']);
% for ichan=1:2
    for jchan = ichan:NumChan
        xdata = InputData(:,ichan);
        ydata = InputData(:,jchan);
        [Pxy(:,ichan,jchan), F] = ...
            cpsd(xdata,ydata,window,NOverlap,FFTLength,SampRate);
        [Cxy,F] = mscohere(xdata,ydata,window,NOverlap,FFTLength,SampRate);
        COH(:,ichan,jchan)=Cxy;
        
        if ichan~=jchan
            Pxy(:,jchan,ichan)=conj(Pxy(:,ichan,jchan));
            COH(:,jchan,ichan)=Cxy;
        end
       
    end
end
return
end
