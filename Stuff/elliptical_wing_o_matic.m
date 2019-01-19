% clc
clear

num_panels = 10;
n = 1;
m = 2;


span = 2.7489;
left = [0 0];
right = [0 span];
rchord = 1;
tchord = 0.01;
% y = linspace(left(2), right(2), num_panels);
y = (((1 - (cos(linspace(pi/2,0,num_panels)))))-1).*-span;
chord = sqrt( (rchord.^2).*(1 - ((y.^2)./(span.^2))));
chord(end) = tchord;

x = (rchord + left(1)) - chord;

y = [-flip(y(2:end)) y]+span;
x = [flip(x(2:end)) x];
chord = [flip(chord(2:end)) chord];

str = sprintf('OPERA INPUT FILE\nALL UNITS IN SI\n\nANALYSIS PARAMETERS\n\nAnalysis Type:	 	strATYPE = 3D\nLifting Surface:	strA2TYPE = THIN\nChordwise Spacing:	strSPACING = NORMAL\n\nMaximum Timesteps: 		valMAXTIME = 0\nTimestep Size: 			valDELTIME = 1\nRelaxed or Fixed Wake: 	strRELAX = FIXED\n\nAngles of Attack: 		seqALPHA = 4\nAngles of Sideslip: 	seqBETA = 0\n\nGEOMETRY\n\nNo. of panels:	valPANELS =	%d\n\nDefines leading edge of wing, all measured in metres:\nKeep vecM the same for all panels on a wing.',2.*(num_panels-1));
for i = 1:(num_panels-1).*2
    panels{i} = sprintf('\nPanel #:%d.\nNumber of spanwise elements:	vecN 		= %d.\nNumber of chordwise elements: 	vecM 		= %d.\nSymmetry about YZ Plane:		strSYM 		= NO\nTrailing edge:					strTE 		= YES\nLeading edge:					strLE       = YES\nxleft		yleft		zleft		chord		epsilon	 	airfoil\n', i, n, m);
    goem = sprintf('%0.5f\t%0.5f\t0\t%0.5f\t0\tNACA 0015\n%0.5f\t%0.5f\t0\t%0.5f\t0\tNACA 0015\n\n',x(i),y(i),chord(i),x(i+1),y(i+1),chord(i+1));
    panels{i} = [panels{i} goem];
    
end

stuff = strjoin(convertCharsToStrings([str panels]));

fid = fopen('G:\GIT\opera\inputs\ellipse.dat','wt');
fprintf(fid, stuff);
fclose(fid);