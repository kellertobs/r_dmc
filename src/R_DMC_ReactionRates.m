%**************************************************************************
%*****  R_DMC REACTION RATES CALCULATION ROUTINE  *************************
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
% FILENAME  R_DMC_ReactionRates.m
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
% HELP:  [Gamma]  =  R_DMC_ReactionRates(VARdis,VAReql,PAR)
% INPUT/OUTPUT: 
% - Input variables structure VARdis contains nt-by-1 row vectors of 
%   disequilibrated model state: f (melt fraction [wt]), and nt-by-nc 
%   matrices Cs (solid composition [wt %]), and Cl (liquid composition [wt %]).
% - Input variables structure VAReql contains nt-by-1 row vectors of 
%   equilibrium model state: f (melt fraction [wt]), and nt-by-nc matrices Cs
%   (solid composition [wt %]), and Cl (liquid composition [wt %]).
% - Input parameters structure PAR contains 1-by-nc row vectors: T0 (melting point
%   [deg C]), A (linear P-coefficient [K/GPa]), B (quadratic P-coefficient
%   [K/GPa^2]), L (latent heat [J/kg]), r (tuning parameter [J/K/kg])
% - Output is 1-by-nc vector Gamma (component melting rates [kg/m3/yr]

function  [Gamma]  =  R_DMC_ReactionRates(VARdis,VAReql,PAR)

%***  get reaction rate factor
R  =  PAR.rho./PAR.tau_react;

%***  get net phase change rate Gamma
GammaNet  =  R .* (VAReql.f - VARdis.f);

%***  get fractional phase compositions
CsFract  =  VARdis.Cl.*PAR.K; CsFract = CsFract./repmat(sum(CsFract,2),1,PAR.nc);
ClFract  =  VARdis.Cs./PAR.K; ClFract = ClFract./repmat(sum(ClFract,2),1,PAR.nc);

%***  set reactive composition to fractional melting/crystallisation
CGamma  =  zeros(size(GammaNet));
%     fractional crystallisation
ind          =  find(repmat(GammaNet,1,PAR.nc) <  0); 
CGamma(ind)  =  CsFract(ind);
%     fractional melting
ind          =  find(repmat(GammaNet,1,PAR.nc) >= 0);
CGamma(ind)  =  ClFract(ind);
%     renormalize reactive composition to unity
CGamma       =  CGamma./repmat(sum(CGamma,2),1,PAR.nc);

%***  get component exchange reaction rates Delta^i
Delta  =  R.*PAR.chi.* ( VAReql.f.*(VAReql.Cl-CGamma) - ...
                         VARdis.f.*(VARdis.Cl-CGamma) );
       
%***  get component reaction rates Gamma^i
Gamma  =  CGamma.*GammaNet + Delta;

end
