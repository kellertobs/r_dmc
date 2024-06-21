%**************************************************************************
%*****  R_DMC EQUILIBRIUM CALCULATION ROUTINES  ***************************
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
% FILENAME  R_DMC_Equilibrium.m
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
% HELP:  [VAR,PAR]  =  R_DMC_Equilibrium(VAR,PAR,type)
% INPUT/OUTPUT: 
% - Input variables structure VAR contains np-by-1 row vectors: P (pressure [GPa]),
%   T (temperature [deg C]), and np-by-nc matrix C (bulk composition [wt %].
% - Input parameters structure PAR contains 1-by-nc row vectors T0 (melting point
%   [deg C]), A (linear P-coefficient [K/GPa]), B (quadratic P-coefficient
%   [K/GPa^2]), L (latent heat [J/kg]), r (tuning parameter [J/K/kg])
% - Input type = 'K': calculate component melting points 'Tm' and 
%   partition coefficients 'K' at given PT-conditions. 
%   Output loaded into PAR structure.
% - Input type = 'T': calculate solidus 'Tsol' and liquidus 'Tliq' at given  
%   bulk composition 'C'. 
%   Output loaded into PAR structure.
% - Input type = 'E': calculate equilibrium melt fraction 'f', and phase
%   compositions 'Cs', 'Cl'. 
%   Output loaded into VAR structure.


function  [VAR,PAR]  =  R_DMC_Equilibrium(VAR,PAR,type)

%*****  main routine to compute various quantities depending input 'type'.

if (~isfield(VAR,'f'))
    VAR.f  = zeros(size(VAR.T));
end
if  (~isfield(VAR,'C'))
    VAR.C = zeros(length(VAR.T),PAR.nc);
end
if (~isfield(VAR,'Cl'))
    VAR.Cl = VAR.C;
end
if (~isfield(VAR,'Cs'))
    VAR.Cs = VAR.C;
end

%*****  compute partition coefficients at P,T *****************************

if strcmp(type,'K')
    
    PAR  =  Tm(VAR.P,PAR);
    PAR  =  K(VAR.T,PAR);
    
%*****  compute Tsol, Tliq at P,C^i  **************************************

elseif strcmp(type,'T')
    
    PAR  =  Tsolidus( VAR,PAR);
    PAR  =  Tliquidus(VAR,PAR);
    
%*****  compute equilibrium f, C_s^i, C_l^i fraction at P,T,C^i  **********

elseif strcmp(type,'E')
    
    VAR  =  Equilibrium(VAR,PAR);
    
end

end


function  [PAR]  =  Tsolidus(VAR,PAR)

%*****  subroutine to compute solidus temperature at given bulk composition

%***  exclude invalid compositions
ii  =  sum(VAR.C(:,sum(VAR.C,1)>1),2)<=1;

%***  get P-dependent pure-component melting Temperature
PAR  =  Tm(VAR.P,PAR);

%***  set starting guess for Tsol
Tsol  =  max(min(min(PAR.Tm)),min(max(max(PAR.Tm)),sum(VAR.C.*PAR.Tm,2)));

%***  get T-dependent partition coefficients Ki
PAR  =  K(Tsol,PAR);

%***  get residual for sum(ci_b/Ki) = 1
r  =  sum(VAR.C./PAR.K,2)-1;

rnorm      =  1;     % initialize residual norm for iterations
n          =  0;     % initialize iteration count
rnorm_tol  =  1e-10; % tolerance for Newton residual
its_tol    =  500;   % maximum number of iterations
eps_T      =  5;     % temperature perturbation for finite differencing, degrees

while rnorm > rnorm_tol  % iterate down to full accuracy
    
    %***  compute partition coefficients Ki at T+eps_T
    PAR  =  K(Tsol+eps_T,PAR);
    
    %***  get residual at T+eps_T
    rp   =  sum(VAR.C./PAR.K,2)-1;
    
    %***  compute partition coefficients Ki at T-eps_T
    PAR  =  K(Tsol-eps_T,PAR);
    
    %***  get residual at T-eps_T
    rm   =  sum(VAR.C./PAR.K,2)-1;
    
    %***  finite difference drdT = (r(T+eps_T)-r(T-eps_T))/2/eps_T
    drdT  =  (rp - rm)./2./eps_T;
    
    %***  apply Newton correction to current guess of Tsol
    Tsol(ii)  =  Tsol(ii) - r(ii)./drdT(ii);
    
    %***  compute partition coefficients Ki at Tsol
    PAR  =  K(Tsol,PAR);
    
    %***  compute residual at Tsol
    r  =  sum(VAR.C./PAR.K,2)-1;
    
    %***  compute Newton residual norm
    rnorm  =  norm(r(ii),2)./sqrt(length(r(ii)));
    
    n  =  n+1;   %  update iteration count
    
    if (n==its_tol)
        error(['!!! Newton solver for solidus T has not converged after ',num2str(its_tol),' iterations !!!']);
    end
    
end

PAR.Tsol  =  Tsol;

end


function  [PAR]  =  Tliquidus(VAR,PAR)

%*****  subroutine to compute liquidus temperature at given bulk composition

%***  exclude invalid compositions
ii  =  sum(VAR.C(:,sum(VAR.C,1)>1),2)<=1;

%***  get P-dependent pure component melting T
PAR  =  Tm(VAR.P,PAR);

%***  set starting guess for Tliq
Tliq  =  max(min(min(PAR.Tm)),min(max(max(PAR.Tm)),sum(VAR.C.*PAR.Tm,2)));

%***  get T-dependent partition coefficients Ki
PAR  =  K(Tliq,PAR);

%***  get residual for sum(ci_b*Ki) = 1
r  =  sum(VAR.C.*PAR.K,2)-1;

rnorm      =  1;     % initialize residual norm for iterations
n          =  0;     % initialize iteration count
rnorm_tol  =  1e-10; % tolerance for Newton residual
its_tol    =  500;   % maximum number of iterations
eps_T      =  1;     % temperature perturbation for finite differencing, degrees

while rnorm > rnorm_tol  % iterate down to full accuracy
    
    %***  compute partition coefficients Ki at T+eps
    PAR  =  K(Tliq+eps_T,PAR);
    
    %***  get residual at T+eps_T
    rp  =  sum(VAR.C.*PAR.K,2)-1;
    
    %***  compute partition coefficients Ki at T-eps_T
    PAR  =  K(Tliq-eps_T,PAR);
    
    %***  get residual at T-eps_T
    rm  =  sum(VAR.C.*PAR.K,2)-1;
    
    %***  compute difference drdT = (r(T+eps_T)-r(T-eps_T))/2/eps_T
    drdT  =  (rp - rm)./2./eps_T;
    
    %***  apply Newton correction to Tliq
    Tliq(ii)  =  Tliq(ii) - r(ii)./drdT(ii);
    
    %***  compute partition coefficients Ki at Tliq
    PAR  =  K(Tliq,PAR);
    
    %***  compute residual at Tliq
    r  =  sum(VAR.C.*PAR.K,2)-1;
    
    %***  compute Newton residual norm
    rnorm  =  norm(r(ii),2)./sqrt(length(r(ii)));
    
    n  =  n+1;   %  update iteration count
    
    if (n==its_tol)
        error(['!!! Newton solver for liquidus T has not converged after ',num2str(its_tol),' iterations !!!']);
    end
    
end

PAR.Tliq  =  Tliq;

end


function  [VAR]  =  Equilibrium(VAR,PAR)

%*****  subroutine to compute equilibrium melt fraction and phase 
%       compositions at given bulk composition, pressure and temperature

%***  get Tsol, Tliq at C = [C1,C2,C3]
PAR  =  Tsolidus( VAR,PAR);
PAR  =  Tliquidus(VAR,PAR);

%***  get P-dependent pure component melting T
PAR  =  Tm(VAR.P,PAR);

%***  get T-dependent partition coefficients Ki
PAR  =  K(max(PAR.Tsol,min(PAR.Tliq,VAR.T)),PAR);

%***  compute residual of unity sum of components
ff   =  repmat(VAR.f,1,PAR.nc);
r    =  sum(VAR.C./(ff+(1-ff).*PAR.K),2) - sum(VAR.C./(ff./PAR.K+(1-ff)),2);

rnorm  =  1;       % initialize residual norm for iterations
n      =  0;       % initialize iteration count
rnorm_tol = 1e-10; % tolerance for Newton residual
its_tol   = 1e3;   % maximum number of iterations

if VAR.T <= PAR.Tsol
    VAR.f(:)  =  0;             % no melt below solidus
elseif VAR.T >= PAR.Tliq
    VAR.f(:)  =  1;             % fully molten above liquidus
else
    VAR.f = max(0,min(1,(VAR.T-PAR.Tsol)./(PAR.Tliq-PAR.Tsol)));
    while rnorm > rnorm_tol     % Newton iteration
        
        %***  compute analytic derivative of residual dr/df
        dr_df  =  - sum(VAR.C.*(1 -PAR.K  )./(ff+(1-ff).*PAR.K).^2,2) ...
                  + sum(VAR.C.*(1./PAR.K-1)./(ff./PAR.K+(1-ff)).^2,2);
        
        %***  apply Newton correction to f
        a = 0.5;
%         while (any((VAR.f - a.*r./dr_df) < -1e-16) || any(VAR.f - a.*r./dr_df > 1-1e-16))
%             a = a/2;
%             if (a<1e-6)
%                 error('R_DMC module: Newton solver for equilibrium melt fraction did not converge with step size %4.4f!',2*a);
%             end
%         end
                
        VAR.f = max(0,min(1,VAR.f - a.*r./dr_df));
                
        
        %***  compute residual of unity sum of components
        ff     =  repmat(VAR.f,1,PAR.nc);
        r      =  sum(VAR.C./(ff+(1-ff).*PAR.K),2) - sum(VAR.C./(ff./PAR.K+(1-ff)),2);
        
        %***  get non-linear residual norm
        rnorm  =  norm(r,2)./sqrt(length(r));
        
        n  =  n+1;  % update iteration count
        if (n==its_tol)
            error(['!!! Newton solver for equilibrium f did not converge after ',num2str(its_tol),' iterations !!!']);
        end
    end
end

%***  safeguard bounds on melt fraction against numerical leaks
VAR.f  =  max(0,min(1,VAR.f));
    
%***  get C_s^i, C_l^i as functions of K^i, C^i and f, safeguard bounds
VAR.Cl  =  max(0,min(1, VAR.C./(ff + (1-ff).*PAR.K) ));
VAR.Cs  =  max(0,min(1, VAR.C./(ff./PAR.K + (1-ff)) ));

end


function  [PAR]  =  Tm(P,PAR)

%***  compute P-dependence of component melting points

switch PAR.Tm_P_mode
    
    case 'quadratic'
        
        %***  Parameterization as in Katz (2003)
        PAR.Tm  =  zeros(size(P,1),PAR.nc);
        for i = 1:PAR.nc
            PAR.Tm(:,i)  =  PAR.T0(i) + PAR.A(i).*P + PAR.B(i).*P.^2;
        end
        
        %***  safeguard: continue melting point with linear slope above Pmax
        Pmax    =  6;
        for i = 1:PAR.nc
            ind = P > Pmax;
            T0   =  PAR.T0(i) + PAR.A(i).*Pmax + PAR.B(i).*Pmax.^2;
            dTdP =  ((PAR.A(i).*Pmax + PAR.B(i).*Pmax.^2)-(PAR.A(i).*(Pmax-0.01) + PAR.B(i).*(Pmax-0.01).^2))./0.01;
            PAR.Tm(ind,i) = T0 + dTdP.*(P(ind)-Pmax);
        end
        
    case 'simonslaw'
        
        %***  Parameterization as in Rudge etal (2011)
        PAR.Tm  =  zeros(size(P,1),PAR.nc);
        for i = 1:PAR.nc
            PAR.Tm(:,i)  =  PAR.T0(i) .* (1 + P/PAR.A(i)) .^ (1/PAR.B(i));
        end
        
end % switch
end % function


function  [PAR]  =  K(T,PAR)

%***  compute T,P-dependence of equilibrium partition coefficients

PAR.L = (T+273.15).*PAR.dS;

switch PAR.K_T_mode
        
    case 'inverse_exp'
        
        %     Parameterization after Rudge, Bercovici, & Spiegelman (2010)
        PAR.K  =  zeros(size(T,1),PAR.nc);
        for i = 1:PAR.nc
            PAR.K(:,i)  =  exp(PAR.L(:,i)./PAR.r(i).*(1./(T+273.15) - 1./(PAR.Tm(:,i)+273.15)));
        end
        
    case 'linear_exp'
        
        %     Linearised exponential dependence
        PAR.K  =  zeros(size(T,1),PAR.nc);
        for i = 1:PAR.nc
            PAR.K(:,i)  =  exp(-PAR.r(i).*(T - PAR.Tm(:,i)));
        end
        
end % switch
end % function
