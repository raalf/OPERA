clear
clc
tic

% Analysis Type and Geometry File

ATYPE = 'LS'; % Lifting Surface
STL = 'CAD Geom/simple_liftingsurface.stl';
% STL = 'CAD Geom/quad.stl';

% STL = 'Cad Geom/lifting_split.stl';

% ATYPE = 'PC'; % Panel Code
% STL = 'CAD Geom/cube.stl';

%% Triangulating Geometry

[TR, ADJE, ELST, VLST, DVE, NELE, EATT, EIDX, ELOC, DNORM, PLEX, DVECT, ALIGN, VATT] = fcnTRIANG(STL, ATYPE);

%% D-Matrix Creation

[D, R] = fcnDWING(EATT, PLEX, NELE, ELOC, ALIGN, VATT);

cellsza = cell2mat(cellfun(@size,VATT,'uni',false));
idxa = cellsza(:,2);
VATT2 = zeros(length(idxa),max(idxa))

%% Plot

[hFig1] = fcnPLOTBODY(1, DVE, NELE, VLST, ELST, DVECT);

%% End

clearvars -except TR ADJE ELST VLST DVE NELE EATT EIDX ELOC DNORM PLEX DVECT ALIGN D R VATT

toc

whos
