%**************************************************************************
%*****  R_DMC USER CONTROL ROUTINE EXAMPLE 1 ******************************
%**************************************************************************
%
% COPIRIGHT (C) 2015
% 	  Tobias Keller [tobias.keller@earth.ox.ac.uk] 
% 	  University of Oxford
% 	  FOALAB
% 	  Earth Sciences
% 	  South Parks Road
% 	  Oxford, OX1 3AN
% 	  United Kingdom
%
% PROJECT   R_DMC
% FILENAME  R_DMC_USER_EX1.m
% 
% LICENSE:  R_DMC is free software: you can redistribute it and/or modify 
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation, version 3.
% 
% R_DMC is distributed in the hope that it will be useful, 
% but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
% GNU General Public License for more details.
% 
% You should have received a copy of the GNU General Public License
% along with this program.  If not, see <http://www.gnu.org/licenses/>.

clear variables;
addpath '../src/' '../cal/'

%**************************************************************************
%*****  USE R_DMC PHASE DIAGRAMS  *****************************************
%**************************************************************************

%*****  SET CALIBRATION FILE FOR R_DMC METHOD  %***************************

R_DMC_CAL_MANTLE

%*****  USER SPECIFIED INPUT PARAMETERS  **********************************

%***  set plotting options   
PAR.np      =  200;    % number of grid points
PLOT.hold   =  false;  % false = clear previous figure; true = hold previous figure
PLOT.LS     =  '-';    % line style
PLOT.scale  =  'lin';  % choose 'log'/'lin' scaling for melt & component fractions

%***  set pressure range [GPa] for phase diagrams  
PAR.Pmin  =  0;
PAR.Pmax  =  5;
PAR.Pref  =  0;

%***  set temperature range [deg C] for phase diagrams  
PAR.Tmin  =  900;
PAR.Tmax  =  1900;
PAR.Tref  =  1350;

%***  set reference composition [wt %] for phase diagrams  
%     provide a ref value for each component
%     components will be normalized to 100%
% PAR.Cref    =  [0,8,16,8,2];
PAR.Cref    =  [0,8,15,10,0.2];
PAR.Cref(1) =  100-sum(PAR.Cref);

%***  set parameters for adiabatic T gradient  
PAR.alpha  =  2e-5;  % thermal expansivity [1/K]
PAR.cp     =  1000;  % specific heat capacity [J/kg/K]
PAR.rho    =  3200;  % ref. rock density [kg/m^3]
PAR.g      =  9.81;  % gravity [m/s^2]

%***  specify phase diagrams to plot  
%     SL-plots: leave 0 for no plots, set to 1 for each required plot
%     others:   leave empty {} for no plots, add row vectors with component 
%               numbers for each required plot 
%               e.g. {[1,2],[2,3],[3,4]}, or {[1,2,3,4],[3,4]}
%     note: only binary and ternary TC-diagrams available
PLOT.SL_plot      =   1;                   % solidus & liquidus: T_sol, T_liq vs P
PLOT.TP_plot      =  {[1,2,3,4,5]};        % melting point: T_m^i vs P
PLOT.KP_plot      =  {[1,2,3,4,5]};        % partition coeff: K^i vs P
PLOT.KT_plot      =  {};        % distriubtion coeff: K^i vs T
PLOT.TC_bin_plot  =  {[4,5]};  % binary  phase loops: T_sol, T_liq vs C^i (binary)
PLOT.TC_ter_plot  =  {};%[1,2,3],[1,3,4],[3,4,5]};    % ternary phase loops: T_sol, T_liq vs C^i (ternary)

%***  specify melting and compositional evolution diagrams to plot  
%     FT- and FP-plots: leave 0 for no plots, set to 1 for each required plot
%     CT- and CP-plots: leave empty {} for no plots, add row vectors with 
%                       component numbers for each required plot 
%                       e.g. {[1,2],[2,3],[3,4]}, or {[1,2,3,4],[3,4]}
PLOT.FT_plot  =   0;             % isobaric  melting: f vs T
PLOT.FP_plot  =   1;             % adiabatic melting: f vs P
PLOT.CT_plot  =  {};%{[1,2,3,4,5]};  % isobaric  composition: C_l^i, C_s^i vs T
PLOT.CP_plot  =  {[1,2,3,4,5]};  % adiabatic composition: C_l^i, C_s^i vs P


%*****  CALL R_DMC PHASE DIAGRAM ROUTINE  *********************************

[VAR,PAR]  =  R_DMC_PhaseDiagrams(PAR,PLOT);

%**************************************************************************

