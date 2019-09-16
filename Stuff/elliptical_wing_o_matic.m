function [] = elliptical_wing_o_matic(num_panels, vecM, valALPHA, strSPANSPACING, strCHORDSPACING, valDELTIME, valMAXTIME, AR, xtcr, vecSYM)

vecN = 1;
rchord = 1;
tchord = 0.1;

syms b
A = ((pi.*rchord.*(0.5.*b))/2);

hspan = double(solve((b.^2)/A == AR, b))./2;

left = [0 0];
right = [0 hspan];

valSPAN = hspan.*2;


if strcmpi(strSPANSPACING,'HALFCOSINE')
    y = (((1 - (cos(linspace(pi/2,0,num_panels)))))-1).*-hspan;
else
    y = linspace(left(2), right(2), num_panels);
end

chord = sqrt( (rchord.^2).*(1 - ((y.^2)./(hspan.^2))));
chord(end) = tchord;

x = (rchord + left(1)) - chord.*xtcr;

y = [-flip(y(2:end)) y];
x = [flip(x(2:end)) x];
chord = [flip(chord(2:end)) chord];

% valAREA = (pi*hspan*rchord)./2;
valAREA = sum(((chord(1:end-1) + chord(2:end))./2).*diff(y));
if strcmpi(vecSYM, 'YES')
    tpanels = num_panels - 1;
else
    tpanels =  2.*(num_panels - 1);
end

str = sprintf('OPERA INPUT FILE\nALL UNITS IN SI\n\nANALYSIS PARAMETERS\n\nAnalysis Type:	 	strATYPE = WING\nLifting Surface:	strA2TYPE = THIN\nSteady or Unsteady:	strA3TYPE = STEADY\n\nChordwise Spacing:	strSPACING = %s\n\nMaximum Timesteps: 		valMAXTIME = %d\nTimestep Size: 			valDELTIME = %0.2f\nRelaxed or Fixed Wake: 	strRELAX = FIXED\n\nAngles of Attack: 		seqALPHA = %0.2f\nAngles of Sideslip: 	seqBETA = 0\n\nDensity:			valDENSITY = 1\nReference Area:		valAREA = %0.5f\nReference Span:		valSPAN = %0.5f\nGust Mode:			valGUSTMODE = 0\n\nGEOMETRY\n\nNo. of panels:	valPANELS =	%d\n\nDefines leading edge of wing, all measured in metres:\nKeep vecM the same for all panels on a wing.', strCHORDSPACING, valMAXTIME, valDELTIME, valALPHA, valAREA, valSPAN, tpanels);
if strcmpi(vecSYM, 'YES')
    for i = num_panels:(num_panels-1).*2
        panels{i} = sprintf('\nPanel #:%d.\nNumber of spanwise elements:	vecN 		= %d.\nNumber of chordwise elements: 	vecM 		= %d.\nSymmetry about YZ Plane:		strSYM 		= %s\nTrailing edge:					strTE 		= YES\nLeading edge:					strLE       = YES\nSpanwise Element Spacing:		strPSPACE    = NORMAL\nxleft		yleft		zleft		chord		epsilon	 	airfoil\n', i, vecN, vecM, vecSYM);
        goem = sprintf('%0.5f\t%0.5f\t0\t%0.5f\t0\tNACA 0015\n%0.5f\t%0.5f\t0\t%0.5f\t0\tNACA 0015\n\n',x(i),y(i),chord(i),x(i+1),y(i+1),chord(i+1));
        panels{i} = [panels{i} goem];
        
    end
else
    
    for i = 1:(2.*(num_panels - 1))
        panels{i} = sprintf('\nPanel #:%d.\nNumber of spanwise elements:	vecN 		= %d.\nNumber of chordwise elements: 	vecM 		= %d.\nSymmetry about YZ Plane:		strSYM 		= %s\nTrailing edge:					strTE 		= YES\nLeading edge:					strLE       = YES\nSpanwise Element Spacing:		strPSPACE    = NORMAL\nxleft		yleft		zleft		chord		epsilon	 	airfoil\n', i, vecN, vecM, vecSYM);
        goem = sprintf('%0.5f\t%0.5f\t0\t%0.5f\t0\tNACA 0015\n%0.5f\t%0.5f\t0\t%0.5f\t0\tNACA 0015\n\n',x(i),y(i),chord(i),x(i+1),y(i+1),chord(i+1));
        panels{i} = [panels{i} goem];
    end
end
panels = panels(~cellfun('isempty',panels));
stuff = strjoin(convertCharsToStrings([str panels]));

fid = fopen('inputs\ellipse.dat','w');
fprintf(fid, stuff);
fclose(fid);

end