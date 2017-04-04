
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
%* Include-File: 'com_h02.inc'                                         *
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
%*****       LIFTING_LINE - COMMON - BLOCK - DATEI 'COM_H02'       *****
%***                   Elementarfluegel-Geometrien                   ***
%***** ==> NACH JEDER AENDERUNG PROGRAMM VOELLIG NEU UEBERSETZEN ! *****
% *********************************************************************
%   *****************************************************************
%
%***********************************************************************
%* LIFTING_LINE Sub-file 'com_h02.inc'                                 *
%*---------------------------------------------------------------------*
%*                                                                     *
%* Created:                                                            *
%*   09/1986 by Karl Heinz Horstmann, EA, DFVLR, Braunschweig          *
%*                                                                     *
%* Version V2.0, November 2005                                         *
%*---------------------------------------------------------------------*
%* -> Please insert your name (and Company), the date and the type of  *
%*    your modification beneath this line                              *
%*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*
%* 05/2004, Carsten Liersch, AS, DLR Braunschweig                      *
%*           - Set-Up of Version 1.0                                   *
%* 11/2005, Carsten Liersch, AS, DLR Braunschweig                      *
%*           - Included reference area and mean chord of panel         *
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
 v_=whos('h02r_1'); if isempty(v_); global h02r_1; if isempty(h02r_1), h02r_1=zeros(1,P42); end; else; if ~v_.persistent; h02r_1_orig=h02r_1; clear h02r_1; global h02r_1; h02r_1=h02r_1_orig; clear h02r_1_orig; end; end;
 v_=whos('h02r_2'); if isempty(v_); global h02r_2; if isempty(h02r_2), h02r_2=zeros(1,P42); end; else; if ~v_.persistent; h02r_2_orig=h02r_2; clear h02r_2; global h02r_2; h02r_2=h02r_2_orig; clear h02r_2_orig; end; end;
 v_=whos('h02r_3'); if isempty(v_); global h02r_3; if isempty(h02r_3), h02r_3=zeros(1,P42); end; else; if ~v_.persistent; h02r_3_orig=h02r_3; clear h02r_3; global h02r_3; h02r_3=h02r_3_orig; clear h02r_3_orig; end; end;
 v_=whos('h02r_4'); if isempty(v_); global h02r_4; if isempty(h02r_4), h02r_4=zeros(1,P42); end; else; if ~v_.persistent; h02r_4_orig=h02r_4; clear h02r_4; global h02r_4; h02r_4=h02r_4_orig; clear h02r_4_orig; end; end;
 v_=whos('h02r_5'); if isempty(v_); global h02r_5; if isempty(h02r_5), h02r_5=zeros(1,P42); end; else; if ~v_.persistent; h02r_5_orig=h02r_5; clear h02r_5; global h02r_5; h02r_5=h02r_5_orig; clear h02r_5_orig; end; end;
 v_=whos('h02r_6'); if isempty(v_); global h02r_6; if isempty(h02r_6), h02r_6=zeros(1,P42); end; else; if ~v_.persistent; h02r_6_orig=h02r_6; clear h02r_6; global h02r_6; h02r_6=h02r_6_orig; clear h02r_6_orig; end; end;
 v_=whos('h02r_7'); if isempty(v_); global h02r_7; if isempty(h02r_7), h02r_7=zeros(1,P42); end; else; if ~v_.persistent; h02r_7_orig=h02r_7; clear h02r_7; global h02r_7; h02r_7=h02r_7_orig; clear h02r_7_orig; end; end;
 v_=whos('h02r_8'); if isempty(v_); global h02r_8; if isempty(h02r_8), h02r_8=zeros(1,P42); end; else; if ~v_.persistent; h02r_8_orig=h02r_8; clear h02r_8; global h02r_8; h02r_8=h02r_8_orig; clear h02r_8_orig; end; end;
 v_=whos('h02r_9'); if isempty(v_); global h02r_9; if isempty(h02r_9), h02r_9=zeros(1,P42); end; else; if ~v_.persistent; h02r_9_orig=h02r_9; clear h02r_9; global h02r_9; h02r_9=h02r_9_orig; clear h02r_9_orig; end; end;
 v_=whos('h02r_10'); if isempty(v_); global h02r_10; if isempty(h02r_10), h02r_10=zeros(1,P42); end; else; if ~v_.persistent; h02r_10_orig=h02r_10; clear h02r_10; global h02r_10; h02r_10=h02r_10_orig; clear h02r_10_orig; end; end;
%
%
%***********************************
%***** Bedeutung der Variablen *****
%***********************************
%
%  KSI1  :\
%  ETA1  : > Koordinaten der l/4-Linie am rechten Rand des
%  ZETA1 :/  Elementarfluegels
%
%  KSI2  :\
%  ETA2  : > Koordinaten der l/4-Linie am linken Rand des
%  ZETA2 :/  Elementarfluegels
%
%  EL1   :\ Fluegeltiefen am rechten bzw. linken Rand des
%  EL2   :/ Elementarfluegels
%
%  FEF   : Flaeche des Elementarfluegels
%
%  LMUEEF: Bezugstiefe des Elementarfluegels
%
%
%************************************
%***** Aufbau des COMMON-Blocks *****
%************************************
%
% common ::;
%% common/H02R/ KSI1, ETA1,   ZETA1, EL1, KSI2, ETA2, ZETA2, EL2,FEF,  LMUEEF;
%% common/H02R/ h02r_1, h02r_2,   h02r_3, h02r_4, h02r_5, h02r_6, h02r_7, h02r_8,h02r_9,  h02r_10;
%
%
% *********************************************************************
%***************     ENDE DER COMMON - BLOCK - DATEI     ***************
% *********************************************************************
