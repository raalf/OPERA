
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
%* Include-File: 'common.inc'                                          *
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
%*****          LIFTING_LINE - PARAMETER - DATEI 'COMMON'          *****
%***                     Parameter (Feldgroessen)                    ***
%***** ==> NACH JEDER AENDERUNG PROGRAMM VOELLIG NEU UEBERSETZEN ! *****
% *********************************************************************
%   *****************************************************************
%
%***********************************************************************
%* LIFTING_LINE Sub-file 'common.inc'                                  *
%*---------------------------------------------------------------------*
%*                                                                     *
%* Created:                                                            *
%*   09/1986 by Karl Heinz Horstmann, EA, DFVLR, Braunschweig          *
%*                                                                     *
%* Version V2.5, January 2016                                          *
%*---------------------------------------------------------------------*
%* -> Please insert your name (and Company), the date and the type of  *
%*    your modification beneath this line                              *
%*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*
%* 05/2004, Carsten Liersch, AS, DLR Braunschweig                      *
%*           - Set-Up of Version 1.0                                   *
%* 11/2005, Carsten Liersch, AS, DLR Braunschweig                      *
%*           - Changed values of some parameters                       *
%* 03/2007, Carsten Liersch, AS, DLR Braunschweig                      *
%*           - Parameter P11 increased to allow more airfoils on the   *
%*             wing                                                    *
%* 01/2016, Jens Rabe, AS, DLR Braunschweig                            *
%*           - Added parameters (P12-P14) for propeller module         *
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
%
%***********************************
%***** Bedeutung der Parameter *****
%***********************************
%
%  NT_MAX: Maximale Anzahl der beim Programmstart uebergebenen Parameter
%
%  ND_MAX: Maximale Anzahl der Unterverzeichnisse im Dateipfad
%
%  LD_MAX: Maximale Laenge eines Verzeichnisnamens im Dateipfad
%
%  p1    : Maximale Anzahl der Teilfluegel (TF)
%
%  p2    : Maximale Anzahl der zu rechnenden Anstellwinkel
%
%  p3    : Maximale Anzahl der zu rechnenden Einzelpunkte
%
%  p4    : Maximale Anzahl der Elementarfluegel (EF)
%          (Bei Symmetrie: auf einer Seite)
%
%  P42   : Maximale Anzahl der EF auf beiden Seiten
%
%  P43   : Größe der Aufzustellenden Gamma-Matrix
%
%  p6    : Maximale Anzahl der zu erreichenden Soll-CA-Werte
%
%  p7    : Maximale Anzahl der Fluegel
%
%  p8    : Maximale Anzahl der TF in Tiefenrichtung pro Fluegel
%
%  p9    : Maximale Laenge des Dateipfades
%
%  P10   : Maximale Anzahl an Profilen pro Fluegel
%
%  P11   : Maximale Anzahl an verwendeten Schnitten pro Profil
%          (pro Fluegel)
%
%  P12   : Maximale Anzahl an Propellern
%
%  P13   : Maximale Anzahl an Winkeln für Geschwindigkeitsfeld
%
%  P14   : Maximale Anzahl an radialen Positionen für Geschwindigkeitsfeld
%
%
%**********************************
%***** Vorzugebende Parameter *****
%**********************************
%
%
%
%********************************************
%***** Automatisch berechnete Parameter *****
%********************************************
%
%
%
% *********************************************************************
%***************        ENDE DER PARAMETER - DATEI       ***************
% *********************************************************************
