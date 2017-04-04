
%***********************************************************************
%* LIFTING_LINE - A multi-lifting-line method for design and analysis  *
%*                of nonplanar wing-configurations                     *
%*                                                                     *
%* Version V2.5, January 2016                                          *
%*                                                                     *
%* Copyright (C) 1986 - 2016 Karl Heinz Horstmann                      *
%*                                                                     *
%***********************************************************************
%*                                                                     *
%* GENERAL INFORMATION                                                 *
%* -------------------                                                 *
%*                                                                     *
%* Include-File: 'com_h04.inc'                                         *
%*                                                                     *
%* This program is free software; you can redistribute it and/or       *
%* modify it under the terms of the GNU General Public License         *
%* as published by the Free Software Foundation; either version 2      *
%* of the License, or (at your option) any later version.              *
%*                                                                     *
%* This program is distributed in the hope that it will be useful,     *
%* but WITHOUT ANY WARRANTY; without even the implied warranty of      *
%* MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the       *
%* GNU General Public License for more details.                        *
%*                                                                     *
%* You should have received a copy of the GNU General Public License   *
%* along with this program; if not, write to the Free Software         *
%* Foundation, Inc., 59 Temple Place - Suite 330, Boston,              *
%* MA  02111-1307, USA.                                                *
%*                                                                     *
%* In case of questions or errors, please contact:                     *
%* Carsten Liersch                                                     *
%*  German Aerospace Center (DLR)                                      *
%*  Institute of Aerodynamics and Flow Technology                      *
%*  Transport Aircraft                                                 *
%*  Lilienthalplatz 7                                                  *
%*  d - 38108 Braunschweig / Germany                                   *
%*  Tel.  : +49 (0)531 / 295 - 2434                                    *
%*  Fax . : +49 (0)531 / 295 - 2320                                    *
%*  E-Mail: carsten.liersch@dlr.de                                     *
%*                                                                     *
%***********************************************************************

%   *****************************************************************
% *********************************************************************
%*****       LIFTING_LINE - COMMON - BLOCK - DATEI 'COM_H04'       *****
%***                       Gamma-Koeffizienten                       ***
%***** ==> NACH JEDER AENDERUNG PROGRAMM VOELLIG NEU UEBERSETZEN ! *****
% *********************************************************************
%   *****************************************************************
%
%***********************************************************************
%* LIFTING_LINE Sub-file 'com_h04.inc'                                 *
%*---------------------------------------------------------------------*
%*                                                                     *
%* Created:                                                            *
%*   09/1986 by Karl Heinz Horstmann, EA, DFVLR, Braunschweig          *
%*                                                                     *
%* Version V2.2, April 2007                                            *
%*---------------------------------------------------------------------*
%* -> Please insert your name (and Company), the date and the type of  *
%*    your modification beneath this line                              *
%*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*
%* 05/2004, Carsten Liersch, AS, DLR Braunschweig                      *
%*           - Set-Up of Version 1.0                                   *
%* 11/2005, Carsten Liersch, AS, DLR Braunschweig                      *
%*           - Removed coefficients                                    *
%*           - Changed name of COMMON-Block                            *
%* 03/2007, Carsten Liersch, AS, DLR Braunschweig                      *
%*           - Changed name of COMMON-Block                            *
%*                                                                     *
%* mm/yyyy, name, company                                              *
%*          description of changes                                     *
%*                                                                     *
%***********************************************************************
%
%
%***********************
%***** Deklaration *****
%***********************
%
 v_=whos('h04r_1'); if isempty(v_); global h04r_1; if isempty(h04r_1), h04r_1=zeros(1,P42); end; else; if ~v_.persistent; h04r_1_orig=h04r_1; clear h04r_1; global h04r_1; h04r_1=h04r_1_orig; clear h04r_1_orig; end; end;
 v_=whos('h04r_2'); if isempty(v_); global h04r_2; if isempty(h04r_2), h04r_2=zeros(1,P42); end; else; if ~v_.persistent; h04r_2_orig=h04r_2; clear h04r_2; global h04r_2; h04r_2=h04r_2_orig; clear h04r_2_orig; end; end;
 v_=whos('h04r_3'); if isempty(v_); global h04r_3; if isempty(h04r_3), h04r_3=zeros(1,P42); end; else; if ~v_.persistent; h04r_3_orig=h04r_3; clear h04r_3; global h04r_3; h04r_3=h04r_3_orig; clear h04r_3_orig; end; end;
%
%
%***********************************
%***** Bedeutung der Variablen *****
%***********************************
%
%  az  : Faktor des konstanten Gliedes der quadratischen Gleichung fuer
%        die  Wirbelverteilung auf jedem Elementarfluegel
%
%  bz  : Faktor des linearen Gliedes der quadratischen Gleichung fuer
%        die  Wirbelverteilung auf jedem Elementarfluegel
%
%  cz  : Faktor des quadratischen Gliedes der quadratischen Gleichung
%        fuer die Wirbelverteilung auf jedem Elementarfluegel
%
%
%************************************
%***** Aufbau des COMMON-Blocks *****
%************************************
%
% common ::;
%% common/H04R/ az, bz, cz;
%% common/H04R/ h04r_1, h04r_2, h04r_3;
%
%
% *********************************************************************
%***************     ENDE DER COMMON - BLOCK - DATEI     ***************
% *********************************************************************
