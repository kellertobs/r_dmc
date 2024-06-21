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
% HELP:  This file contains a test calibration for idealised, linearised 
% mantle melting in a five-component compositional space, with components 
% c1,...,c5 sorter in order of increasing incompatibility and decreasing 
% concentration.

% set number of components for R_DMC method  [USER SPEC]
PAR.nc  =  5;

% set string with component name to appear in figure legends
PAR.CompStr  =  {'for','fay','opx','cpx','vol'};

% specify calibration parameters for all components below [USER SPEC]

% set pure component melting points T_m^i at P=0
PAR.T0(1)  =  1890;
PAR.T0(2)  =  1205;
PAR.T0(3)  =  1175;
PAR.T0(4)  =  1000;
PAR.T0(5)  =   800;

% choose type of parameterisation for T_m^i(P)
PAR.Tm_P_mode  = 'simonslaw';

% set first coeff. for P-dependence of T_m^i [GPa]
PAR.A(1)   =   10.83;
PAR.A(2)   =   15.78;
PAR.A(3)   =   4.6;
PAR.A(4)   =   1.8;
PAR.A(5)   =   5.6;

% set second coeff. for P-dependence of T_m^i [1]
PAR.B(1)   =  3.7;
PAR.B(2)   =  1.6;
PAR.B(3)   =  2.6;
PAR.B(4)   =  3;
PAR.B(5)   =  2;

% set entropy gain of fusion DeltaS [J/K]
PAR.dS(1)  =  350;
PAR.dS(2)  =  350;
PAR.dS(3)  =  350;
PAR.dS(4)  =  350;
PAR.dS(5)  =  350;

% choose type of parameterisation for K^i(T)
PAR.K_T_mode  = 'inverse_exp';

% set coeff. for T-dependence of partition coefficients K^i [1/K]
PAR.r(1)   =  80;
PAR.r(2)   =  60;
PAR.r(3)   =  8;
PAR.r(4)   =  10;
PAR.r(5)   =  12;


