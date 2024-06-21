%**************************************************************************
%*****  R_DMC CALIBRATION FILE  *******************************************
%**************************************************************************

% set number of components for R_DMC method  [USER SPEC]
PAR.nc  =  4;

% set string with component name to appear in figure legends
PAR.CompStr  =  {'dunite~','morb~','hmorb~','cmorb~'};

% specify calibration parameters for all components below [USER SPEC]

% set pure component melting points T_m^i at P=0
PAR.T0(1)  =  1780;                  
PAR.T0(2)  =  976;
PAR.T0(3)  =  651;
PAR.T0(4)  =  575;

% choose type of parameterisation for T_m^i(P)
PAR.Tm_P_mode  = 'quadratic';

% set coeff. for linear P-dependence of T_m^i [K/GPa]
PAR.A(1)   =   45;
PAR.A(2)   =  107.4;
PAR.A(3)   =   32.2;
PAR.A(4)   =   22.1;

% set coeff. for quadratic P-dependence of T_m^i [K/GPa^2]
PAR.B(1)   =  -2;
PAR.B(2)   =  -2.83;
PAR.B(3)   =  -1.20;
PAR.B(4)   =  -1.42;

% set latent heat of pure components L^i [J/kg]
PAR.L(1)   =  600e3;
PAR.L(2)   =  450e3;
PAR.L(3)   =  350e3;
PAR.L(4)   =  350e3;

% choose type of parameterisation for K^i(T)
PAR.K_T_mode  = 'inverse_exp';

% set coeff. for T-dependence of distribution coefficients K^i
PAR.r(1)   =  70;
PAR.r(2)   =  40;
PAR.r(3)   =  40;
PAR.r(4)   =  40;


