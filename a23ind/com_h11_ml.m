
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
%* Include-File: 'com_h11.inc'                                         *
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
%*****       LIFTING_LINE - COMMON - BLOCK - DATEI 'COM_H11'       *****
%***                  Steuergroessen (Kaneale, Zeit)                 ***
%***** ==> NACH JEDER AENDERUNG PROGRAMM VOELLIG NEU UEBERSETZEN ! *****
% *********************************************************************
%   *****************************************************************
%
%***********************************************************************
%* LIFTING_LINE Sub-file 'com_h11.inc'                                 *
%*---------------------------------------------------------------------*
%*                                                                     *
%* Created:                                                            *
%*   09/1986 by Karl Heinz Horstmann, EA, DFVLR, Braunschweig          *
%*                                                                     *
%* Version V2.3, October 2010                                          *
%*---------------------------------------------------------------------*
%* -> Please insert your name (and Company), the date and the type of  *
%*    your modification beneath this line                              *
%*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*
%* 05/2004, Carsten Liersch, AS, DLR Braunschweig                      *
%*           - Set-Up of Version 1.0                                   *
%* 11/2005, Carsten Liersch, AS, DLR Braunschweig                      *
%*           - Splitted COMMON-Blocks according to types               *
%* 03/2007, Carsten Liersch, AS, DLR Braunschweig                      *
%*           - New variables for global warning handling and           *
%*             command-line option 'Quiet'                             *
%* 10/2010, Carsten Liersch, AS, DLR Braunschweig                      *
%*           - New variable for identification, whether projected area *
%*             and spanwidth are lying in XY- or in XZ-plane           *
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
 v_=whos('h11i_1'); if isempty(v_); global h11i_1; if isempty(h11i_1), h11i_1=0; end; else; if ~v_.persistent; h11i_1_orig=h11i_1; clear h11i_1; global h11i_1; h11i_1=h11i_1_orig; clear h11i_1_orig; end; end;
 v_=whos('h11i_2'); if isempty(v_); global h11i_2; if isempty(h11i_2), h11i_2=0; end; else; if ~v_.persistent; h11i_2_orig=h11i_2; clear h11i_2; global h11i_2; h11i_2=h11i_2_orig; clear h11i_2_orig; end; end;
 v_=whos('h11i_3'); if isempty(v_); global h11i_3; if isempty(h11i_3), h11i_3=0; end; else; if ~v_.persistent; h11i_3_orig=h11i_3; clear h11i_3; global h11i_3; h11i_3=h11i_3_orig; clear h11i_3_orig; end; end;
 v_=whos('h11i_4'); if isempty(v_); global h11i_4; if isempty(h11i_4), h11i_4=0; end; else; if ~v_.persistent; h11i_4_orig=h11i_4; clear h11i_4; global h11i_4; h11i_4=h11i_4_orig; clear h11i_4_orig; end; end;
 v_=whos('h11i_5'); if isempty(v_); global h11i_5; if isempty(h11i_5), h11i_5=0; end; else; if ~v_.persistent; h11i_5_orig=h11i_5; clear h11i_5; global h11i_5; h11i_5=h11i_5_orig; clear h11i_5_orig; end; end;
 v_=whos('h11i_6'); if isempty(v_); global h11i_6; if isempty(h11i_6), h11i_6=0; end; else; if ~v_.persistent; h11i_6_orig=h11i_6; clear h11i_6; global h11i_6; h11i_6=h11i_6_orig; clear h11i_6_orig; end; end;
 v_=whos('h11i_7'); if isempty(v_); global h11i_7; if isempty(h11i_7), h11i_7=0; end; else; if ~v_.persistent; h11i_7_orig=h11i_7; clear h11i_7; global h11i_7; h11i_7=h11i_7_orig; clear h11i_7_orig; end; end;
 v_=whos('h11i_8'); if isempty(v_); global h11i_8; if isempty(h11i_8), h11i_8=0; end; else; if ~v_.persistent; h11i_8_orig=h11i_8; clear h11i_8; global h11i_8; h11i_8=h11i_8_orig; clear h11i_8_orig; end; end;
 v_=whos('h11i_9'); if isempty(v_); global h11i_9; if isempty(h11i_9), h11i_9=0; end; else; if ~v_.persistent; h11i_9_orig=h11i_9; clear h11i_9; global h11i_9; h11i_9=h11i_9_orig; clear h11i_9_orig; end; end;
 v_=whos('h11l_1'); if isempty(v_); global h11l_1; if isempty(h11l_1), h11l_1=false; end; else; if ~v_.persistent; h11l_1_orig=h11l_1; clear h11l_1; global h11l_1; h11l_1=h11l_1_orig; clear h11l_1_orig; end; end;
 v_=whos('h11l_2'); if isempty(v_); global h11l_2; if isempty(h11l_2), h11l_2=false; end; else; if ~v_.persistent; h11l_2_orig=h11l_2; clear h11l_2; global h11l_2; h11l_2=h11l_2_orig; clear h11l_2_orig; end; end;
 v_=whos('h11l_3'); if isempty(v_); global h11l_3; if isempty(h11l_3), h11l_3=false; end; else; if ~v_.persistent; h11l_3_orig=h11l_3; clear h11l_3; global h11l_3; h11l_3=h11l_3_orig; clear h11l_3_orig; end; end;
 v_=whos('h11l_4'); if isempty(v_); global h11l_4; if isempty(h11l_4), h11l_4=false; end; else; if ~v_.persistent; h11l_4_orig=h11l_4; clear h11l_4; global h11l_4; h11l_4=h11l_4_orig; clear h11l_4_orig; end; end;
 v_=whos('h11l_5'); if isempty(v_); global h11l_5; if isempty(h11l_5), h11l_5=false; end; else; if ~v_.persistent; h11l_5_orig=h11l_5; clear h11l_5; global h11l_5; h11l_5=h11l_5_orig; clear h11l_5_orig; end; end;
 v_=whos('h11r_1'); if isempty(v_); global h11r_1; if isempty(h11r_1), h11r_1=0; end; else; if ~v_.persistent; h11r_1_orig=h11r_1; clear h11r_1; global h11r_1; h11r_1=h11r_1_orig; clear h11r_1_orig; end; end;
 v_=whos('h11r_2'); if isempty(v_); global h11r_2; if isempty(h11r_2), h11r_2=0; end; else; if ~v_.persistent; h11r_2_orig=h11r_2; clear h11r_2; global h11r_2; h11r_2=h11r_2_orig; clear h11r_2_orig; end; end;
 v_=whos('h11r_3'); if isempty(v_); global h11r_3; if isempty(h11r_3), h11r_3=zeros(1,3); end; else; if ~v_.persistent; h11r_3_orig=h11r_3; clear h11r_3; global h11r_3; h11r_3=h11r_3_orig; clear h11r_3_orig; end; end;
%
%
%***********************************
%***** Bedeutung der Variablen *****
%***********************************
%
%  IERGEB     : Ausgabedatei-Kanal fuer die Gesamtbeiwerte (TECPLOT)
%
%  IERGEB1    : Ausgabedatei-Kanal fuer Beiwertverteilungen (TECPLOT)
%
%  IERGEB2    : Ausgabedatei-Kanal fuer Zirkulationsverteilungen (TECPLOT)
%
%  IERGEB3    : Ausgabedatei-Kanal fuer die Flugzeuggeometrie (TECPLOT)
%
%  IERGEB4    : Ausgabedatei-Kanal fuer die Einzelpunkte (TECPLOT)
%
%  IDRUCK     : Ausgabedatei-Kanal fuer die allgemeinen Ergebnisse
%
%  IMONITOR   : Ausgabedatei-Kanal fuer Kontrollausgaben
%
%  IXMLOUT    : Ausgabedatei-Kanal fuer die XML-Ausgabe
%
%  DB_LEVEL   : Debug-Level (Hoehere Zahl bedeutet: mehr Ausgaben in die
%               Monitordatei)
%
%  QUIET_MODE : true_ml  = Standard-Bildschirmausgabe unterdruecken
%               false_ml = Normale Standard-Bildschirmausgabe
%
%  TEC_OUT    : true_ml  = TECPLOT-Ausgaben aktivieren
%               false_ml = TECPLOT-Ausgaben deaktivieren
%
%  WARN       : Marker-Flag zur Anzeige, ob waehrend des Programmlaufs
%               bereits Warnungen aufgetreten sind
%               true_ml  = Warnungen aufgetreten
%               false_ml = Keine Warnungen aufgetreten
%
%  XML_OUT    : true_ml  = XML-Ausgaben aktivieren
%               false_ml = XML-Ausgaben deaktivieren
%  PROJ_XY    : Marker-Flag zur Anzeige, ob die projizierte Spannweite
%               und Flaeche in der xy-Ebene oder in der XZ-Ebene liegen
%               true_ml  = Projektion in die XY-Ebene (Normalfall)
%               false_ml = Projektion in die XZ-Ebene (Wenn die Spannweite
%                       in der XY-Projektion zu klein waere)
%
%  A_TIMER    : Zeitlaufvariable - Anfang der Zeitmessstrecke
%
%  E_ETIMER   : Zeitlaufvariable - Ende der Zeitmessstrecke
%
%  TIMER (1)  : Summe aller Zeiten fuer Vorlaufrechnungen
%        (2)  : Summe aller Zeiten fuer Loesung des Systems
%        (3)  : Summe aller Zeiten fuer Nachlaufrechnungen
%
%
%
%************************************
%***** Aufbau des COMMON-Blocks *****
%************************************
%
% common ::;
%% common/H11I/ IERGEB,   IERGEB1, IERGEB2, IERGEB3, IERGEB4, IDRUCK,IMONITOR, IXMLOUT, DB_LEVEL;
%% common/H11I/ h11i_1,   h11i_2, h11i_3, h11i_4, h11i_5, h11i_6,h11i_7, h11i_8, h11i_9;
% common ::;
%% common/H11L/ QUIET_MODE, TEC_OUT, WARN, XML_OUT, PROJ_XY;
%% common/H11L/ h11l_1, h11l_2, h11l_3, h11l_4, h11l_5;
% common ::;
%% common/H11R/ A_TIMER, E_TIMER, TIMER;
%% common/H11R/ h11r_1, h11r_2, h11r_3;
%
% *********************************************************************
%***************     ENDE DER COMMON - BLOCK - DATEI     ***************
% *********************************************************************
