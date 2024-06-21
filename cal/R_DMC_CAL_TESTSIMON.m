%**************************************************************************
%*****  R_DMC CALIBRATION EXAMPLE FILE  ***********************************
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
% FILENAME  R_DMC_Cal_HydrCarbMantle.m
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
% 
%
% HELP:  This file contains a test calibration for hydrated mantle melting 
% in a three-component compositional space of dunite, MORB, and hydrated 
% silicate with 5 wt % H2O.

% set number of components for R_DMC method  [USER SPEC]
PAR.nc  =  3;

% set string with component name to appear in figure legends
PAR.CompStr  =  {'dun~','pxn~','hsi~'};

% specify calibration parameters for all components below [USER SPEC]

% set pure component melting points T_m^i at P=0
PAR.T0(1)  =  1800;                  
PAR.T0(2)  =  1000;
PAR.T0(3)  =  700;

% choose type of parameterisation for T_m^i(P)
PAR.Tm_P_mode  = 'simonslaw';

% set first coeff. for P-dependence of T_m^i [GPa]
PAR.A(1)   =   10;
PAR.A(2)   =   3;
PAR.A(3)   =   8;

% set second coeff. for P-dependence of T_m^i [1]
PAR.B(1)   =  3.8;
PAR.B(2)   =  3.3;
PAR.B(3)   =  3.5;

% set latent heat of pure components L^i [J/kg]
PAR.L(1)   =  600e3;
PAR.L(2)   =  450e3;
PAR.L(3)   =  350e3;

% choose type of parameterisation for K^i(T)
PAR.K_T_mode  = 'inverse_exp';

% set coeff. for T-dependence of distribution coefficients K^i
PAR.r(1)   =  50;
PAR.r(2)   =  30;
PAR.r(3)   =  24;


