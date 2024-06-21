%**************************************************************************
%*****  R_DMC USER CONTROL ROUTINE EXAMPLE 2 ******************************
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
% FILENAME  R_DMC_USER_EX2.m
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
%*****  USE R_DMC REACTIVE EQUILIBRATION  *********************************
%**************************************************************************

%*****  SET CALIBRATION FILE FOR R_DMC METHOD  %***************************

R_DMC_CAL_JPET16_REF

%*****  USER SPECIFIED INPUT PARAMETERS  **********************************

%***  set plotting options   
PAR.nt      =  500;    % number of time steps
PLOT.hold   =  false;  % false = clear previous figure; true = hold previous figure
PLOT.LS     =  '-';    % line style
PLOT.scale  =  'lin';  % choose 'log'/'lin' scaling for melt & component fractions

%***  set reactive model parameters
PAR.tau_react  =  1e4;  % [yr] reaction time constant
PAR.chi        =  1;    % [1] assimilation ratio

%***  set time span for reactive equilibration
PAR.tend  =  1e6;  % [yr]

%***  set initial pressure at t=0 [GPa]
PAR.Pstart  =  4;  % [GPa]

%***  set de/compression rate applied during reactive equilibration
PAR.PRate  =  -4/PAR.tend;  % [GPa/yr]

%***  set initial potential temperature at t=0 
PAR.Tstart  =  1350;  % [deg C]

%***  set heating/cooling rate applied during reactive equilibration
%     this rate is applied in addition to latent and adiabatic heat exchange rates
PAR.TRate  =  0/PAR.tend;  % [deg C/yr]

%***  set initial bulk composition at t=0 
%     provide a ref value for each component
%     components will be normalized to 100%
PAR.Cstart  =  [70,30,0.2,0.05];  % [wt %]

%***  set melt addition/removal rate applied during reactive equilibration 
%     this rate is applied in addition to reaction rate
PAR.fRate  =  0/PAR.tend;  % [wt/yr]

%***  set parameters for adiabatic T gradient  
PAR.alpha  =  2e-5;  % thermal expansivity [1/K]
PAR.cp     =  1000;  % specific heat capacity [J/kg/K]
PAR.rho    =  3200;  % ref. rock density [kg/m^3]
PAR.g      =  9.81;  % gravity [m/s^2]

%***  specify reactive equilibration diagrams to plot
%     FT- and FP-plots: leave 0 for no plots, set to 1 for each required plot
%     CT- and CP-plots: leave empty {} for no plots, add row vectors for each required plot 
%                       e.g. {[1,2],[2,3],[3,4]}, or {[1,2,3,4],[3,4]}
PLOT.PT_plot  =   1;             % PT-path evolution with time: P vs T
PLOT.F_plot   =   1;             % melt evolution: f vs time
PLOT.R_plot   =  {[1,2],[3,4]};  % component reaction rates evolution: Gamma^i vs time
PLOT.C_plot   =  {[1,2],[3,4]};  % evolution of disequilibrium composition: C_s,l^i vs time
PLOT.CB_plot  =  {[1,2],[3,4]};  % evolution of mixture composition: C_b^i vs time
PLOT.EQ_plot  =  1;              % show equilibrium f, C_s,f^i for comparison
PLOT.FR_plot  =  1;              % show fractional C_s,l^i for comparison


%*****  CALL R_DMC REACTIVE EQUILIBRATION ROUTINE  ************************

[VAR,PAR]  =  R_DMC_ReactiveEquilibration(PAR,PLOT);

%**************************************************************************

