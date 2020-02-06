function FFTParam = FFTParameters
%
% a function to create a FFTParam structure for postprocessing of Nicolet data files
%
% use:	FFTParam = FFTParameters
%
%	Inputs:

%	Outputs:	FFTParam - structure containing the following  --structure--
%
%			FFTLength - Length of the FFT
%			SampRate - Data Acquistion Sample Rate - (note:  will eventually come from logdata)
%			WinType - Data analysis window type (ie - hanning)
%			PercOverlap - Data analysis percent overlap
%			DFLAG - Data detrending flag (set to 'none')

%
%	DWDeVilbiss		24-Aug-1999				New
%
% pop a GUI asking for the FFT parameters if not specified
for i=1:100000
    dummy=log(i);
end
prompt={'Enter the FFT Length:',...
      'Enter the Sample Rate:',...
      'Enter the Window Type:',...
      'Enter the Percent Overlap:', ...
      'Enter the Detrend Mode:'};
def={'2048','2048','hann','50','none'};
title='Input for Power Spectral Density Calculations';
lineNo=1;

answer=inputdlg(prompt,title,lineNo,def);
if length(answer)~= 0
   FFTParam.FFTLength = str2num(char(answer{1}));
   FFTParam.SampRate = str2num(char(answer{2}));
   FFTParam.WinType = char(answer{3});
   FFTParam.PercOverlap = str2num(char(answer{4}));
   FFTParam.DFlag = char(answer{5});
end
