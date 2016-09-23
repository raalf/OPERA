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

[TR, ADJE, ELST, VLST, DVE, NELE, EATT, EIDX, ELOC, ...
    PLEX, DVECT, ALIGN, VATT, VNORM, CENTER] = fcnIMPORTGEOM(STL, ATYPE);

%% D-Matrix Creation

[D, R] = fcnDWING(EATT, PLEX, NELE, ELOC, ALIGN, VLST, VNORM, CENTER, DVE, DVECT);

%% Plot

[hFig1] = fcnPLOTBODY(1, DVE, NELE, VLST, ELST, DVECT, CENTER);

%% End

clearvars -except TR ADJE ELST VLST DVE NELE EATT EIDX ELOC DNORM PLEX DVECT ALIGN D R VATT VNORM CENTER

toc

whos
