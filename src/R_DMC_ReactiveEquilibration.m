%**************************************************************************
%*****  R_DMC REACTIVE EQUILIBRATION ROUTINE  *****************************
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
% FILENAME  R_DMC_ReactiveEquilibration.m
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
% HELP:  [VAR,PAR] = R_DMC_ReactiveEquilibration(PAR,PLOT)
% INPUT/OUTPUT: 
% - Input parameters structure PAR:
%   T0      melting point at P=0     [deg C]    1-by-nc row vector
%   A       linear P-coeff           [K/GPa]    "
%   B       quadratic P-coeff        [K/GPa^2]  "
%   L       latent heat              [J/kg]     " 
%   r       tuning parameter         [J/K/kg]   "
%   nc      number of components     [1]        integer
%   nt      number of timesteps      [1]        integer
%   tau     reaction time constant   [yr]       scalar
%   chi     assimilation ratio       [1]        scalar
%   tend    equilibration end time   [yr]       scalar
%   Pstart  initial pressure         [GPa]      scalar
%   PRate   de/compression rate      [GPa/yr]   scalar
%   Tstart  initial temperature      [deg C]    scalar
%   TRate   heating/cooling rate     [deg C]    scalar
%   fRate   melt add./removal rate   [wt/yr]    scalar
%   Cstart  reference bulk comp.     [wt %]     1-by-nc row vector
%   alpha thermal expansivity        [1/K]      scalar
%   cp    specific heat              [J/kg/K]   scalar
%   rho   density                    [kg/m3]    scalar
%   g     gravity                    [m/s2]     scalar
% - Input plotting control structure PLOT:
%   hold     hold previous figures     [true/false]        logical
%   LS       plotting line style       ['-',':',etc]       string
%   scale    choice of plot scaling    ['lin'/'log']       string
%   PT_plot  plot PT-path              [1/0]               integer
%   F_plot   plot f vs t               [1/0]               integer
%   R_plot   plot Gamma^i vs t         e.g. {[1,2],[3,4]}  cell array
%   C_plot   plot C_s,l^i vs t         e.g. {[1,2],[3,4]}  cell array
%   CB_plot  plot bulk C^i vs t        e.g. {[1,2],[3,4]}  cell array
%   CE_plot  plot equil. C_s,l^i vs t  e.g. {[1,2],[3,4]}  cell array
%   CF_plot  plot fract. C_s,l^i vs t  e.g. {[1,2],[3,4]}  cell array
% - Output loaded into VAR, PAR structures.
%
% ODE solver contained in this routine solves the following system of
% coupled equations for temperature T [deg C], pressure [GPa], melt
% fraction [wt], and phase compositions Cs, Cl [wt]. Reactive equilibration
% driven by user specified disequilibration rates for de/compression "PRate",
% heating/cooling "TRate" and de/compaction "fRate".
%
% (1)   dP/dt   =  PRate
%
% (2)   dT/dt   =  TRate + ( alpha.T.PRate - Sum_i(Gamma^i.L^i) ) / (rho.cp)
%
% (3)   df/dt   =  fRate + SUM_i(Gamma^i) / rho
%
% (4a)  dCs/dt  =  -(Gamma^i - Cs^i.SUM_i(Gamma^i)) / ((1-f).rho)
% (4b)  dCl/dt  =   (Gamma^i - Cl^i.SUM_i(Gamma^i)) / (   f .rho)


function  [VAR,PAR] = R_DMC_ReactiveEquilibration(PAR,PLOT)

%*****  set up equilibration by reaction for specified input  *************

%***  initialize vector to store component reaction rates
PAR.Gamma  =  zeros(1,PAR.nc);

%***  normalize input C_ref to 1
PAR.Cstart  =  PAR.Cstart./repmat(sum(PAR.Cstart),1,PAR.nc);

%***  set up starting point model state P,T,c
VAR.P  =  PAR.Pstart;
VAR.T  =  PAR.Tstart.*exp(PAR.Pstart.*1e9.*PAR.alpha./PAR.rho./PAR.cp);
VAR.C  =  PAR.Cstart;

%***  call R_DMC for equilibrium f,C_s^i,C_l^i
[VAR,~]  =  R_DMC_Equilibrium(VAR,PAR,'E');

%***  set initial solution vector Y0 = [f,C_s^i,C_l^i,Gamma^i]
Y0  =  [VAR.P;VAR.T./1000;VAR.f;VAR.Cs(:);VAR.Cl(:);PAR.Gamma(:)];


%*****  compute equilibration by reaction  ********************************

opt           =  odeset('RelTol',1e-9,'NonNegative',1:length(Y0)-PAR.nc);
[PAR.time,Y]  =  ode15s(@RatesOfChange,linspace(0,PAR.tend,PAR.nt),Y0,opt,PAR);

%***  read solution from ODE solver into data structure
VAR.P      =  Y(:,1);  % [GPa]
VAR.T      =  Y(:,2).*1000;  % [deg C]
VAR.f      =  Y(:,3);  % [wt]
VAR.Cs     =  Y(:,3+       (1:PAR.nc));  % [wt]
VAR.Cl     =  Y(:,3+PAR.nc+(1:PAR.nc));  % [wt]
PAR.Gamma  =  Y(:,end-PAR.nc+1:end).*1000;  % [g/m3/yr]

%***  post-process bulk, equilibrium, and fractional phase compositions
VAR.C      =  repmat(VAR.f,1,PAR.nc).*VAR.Cl + (1-repmat(VAR.f,1,PAR.nc)).*VAR.Cs;  % [wt]
[~  ,PAR]  =  R_DMC_Equilibrium(VAR,PAR,'K');
[VAReq,~]  =  R_DMC_Equilibrium(VAR,PAR,'E');
VAR.fEq    =  VAReq.f;  VAR.CsEq   =  VAReq.Cs;  VAR.ClEq   =  VAReq.Cl;  clear VAReq;
VAR.CsFr   =  VAR.Cl.*PAR.K; VAR.CsFr = VAR.CsFr./repmat(sum(VAR.CsFr,2),1,PAR.nc);
VAR.ClFr   =  VAR.Cs./PAR.K; VAR.ClFr = VAR.ClFr./repmat(sum(VAR.ClFr,2),1,PAR.nc);


%*****  plot reactive equilibration results  ******************************

%***  initialize some items for plotting
load R_DMC_Colormap;  % load custom colormap
f        =  1;        % initialize figure counter
fs       =  {'FontSize',22};
ls       =  {'LineStyle',PLOT.LS};
lw       =  {'LineWidth',2};
loc      =  {'Location','north'};
ort      =  {'Orientation','Horizontal'};
tx       =  {'Interpreter','Latex'};
tl       =  {'TickLabelInterpreter','Latex'};
holdfig  =  PLOT.hold;
scaling  =  PLOT.scale;
color    =  [flipud(R_DMC_Colormap(8:16:end,:));flipud(R_DMC_Colormap(12:16:end,:))];
set(0,'DefaultFigureColor',[1,1,1]);
set(0,'DefaultAxesLineWidth',1.5);
set(0,'DefaultFigureRenderer','zbuffer');

%***  plot PT-path with time
if (~isempty(PLOT.PT_plot))
    figure(f); f=f+1;
    if ~holdfig; clf; end; hold on;
    plot(VAR.T,VAR.P,ls{:},lw{:},'Color',color(4,:));
    box on; grid on; axis ij;
    set(gca,fs{:},tl{:});
    xlabel('Temperature [deg C]',fs{:},tx{:});
    ylabel('Pressure [GPa]',fs{:},tx{:});
    drawnow;
    if ~holdfig; hold off; end
end

%***  plot melt evolution with time
if (~isempty(PLOT.F_plot))
    figure(f); f=f+1;
    if ~holdfig; clf; end; hold on;
    if strcmp(scaling,'lin'); plot(PAR.time./1000,               VAR.f  ,ls{:},lw{:},'Color',color(1,:)); end;
    if strcmp(scaling,'log'); plot(PAR.time./1000,log10(max(1e-3,VAR.f)),ls{:},lw{:},'Color',color(1,:)); end;
    if PLOT.EQ_plot
        if strcmp(scaling,'lin'); plot(PAR.time./1000,               VAR.fEq  ,':',lw{1},0.5,'Color',color(1,:)); end;
        if strcmp(scaling,'log'); plot(PAR.time./1000,log10(max(1e-3,VAR.fEq)),':',lw{1},0.5,'Color',color(1,:)); end;
    end
    box on; grid on;
    set(gca,'xlim',[0,PAR.tend./1000],fs{:},tl{:});
    xlabel('Time [kyr]',fs{:},tx{:});
    ylabel('Melt fraction [wt \%]',fs{:},tx{:});
    drawnow;
    if ~holdfig; hold off; end
end

%***  plot component reaction rates with time
for k = 1:size(PLOT.R_plot,2)
    figure(f); f=f+1;
    if ~holdfig; clf; end; hold on;
    for i = PLOT.R_plot{k}(:)'
        if strcmp(scaling,'lin'); plot(PAR.time./1000,      PAR.Gamma(:,i) ,ls{:},lw{:},'Color',color(i,:)); end;
        if strcmp(scaling,'log'); plot(PAR.time./1000,log10(PAR.Gamma(:,i)),ls{:},lw{:},'Color',color(i,:)); end;
    end
    box on; grid on;
    set(gca,'xlim',[0,PAR.tend./1000],fs{:},tl{:});
    leg = legend(PAR.CompStr{PLOT.R_plot{k}(:)'}); set(leg,loc{:},ort{:},fs{:},tx{:});
    xlabel('Time [kyr]',fs{:},tx{:});
    ylabel('Reaction rates [g/m$^3$/yr]',fs{:},tx{:});
    drawnow;
    if ~holdfig; hold off; end
end

%***  plot compositions with time
for k = 1:size(PLOT.C_plot,2)
    % plot solid composition
    figure(f); f=f+1;
    if ~holdfig; clf; end; hold on;
    for i = PLOT.C_plot{k}(:)'
        if strcmp(scaling,'lin'); plot(PAR.time./1000,      VAR.Cs(:,i) ,ls{:},lw{:},'Color',color(i,:)); end;
        if strcmp(scaling,'log'); plot(PAR.time./1000,log10(VAR.Cs(:,i)),ls{:},lw{:},'Color',color(i,:)); end;
    end
    % add equilibrium solid composition if specified
    if PLOT.EQ_plot
        for i = PLOT.R_plot{k}(:)'
            if strcmp(scaling,'lin'); plot(PAR.time./1000,      VAR.CsEq(:,i) ,':',lw{1},0.5,'Color',color(i,:)); end;
            if strcmp(scaling,'log'); plot(PAR.time./1000,log10(VAR.CsEq(:,i)),':',lw{1},0.5,'Color',color(i,:)); end;
        end
    end
    % add fractional solid composition if specified
    if PLOT.FR_plot
        for i = PLOT.R_plot{k}(:)'
            if strcmp(scaling,'lin'); plot(PAR.time./1000,      VAR.CsFr(:,i) ,'-.',lw{1},0.5,'Color',color(i,:)); end;
            if strcmp(scaling,'log'); plot(PAR.time./1000,log10(VAR.CsFr(:,i)),'-.',lw{1},0.5,'Color',color(i,:)); end;
        end
    end
    box on; grid on;
    set(gca,'xlim',[0,PAR.tend./1000],fs{:},tl{:});
    leg = legend(PAR.CompStr{PLOT.C_plot{k}(:)'}); set(leg,loc{:},ort{:},fs{:},tx{:});
    xlabel('Time [kyr]',fs{:},tx{:});
    ylabel('Rock composition [wt \%]',fs{:},tx{:});
    drawnow;
    if ~holdfig; hold off; end
    % plot liquid composition
    figure(f); f=f+1;
    if ~holdfig; clf; end; hold on;
    for i = PLOT.C_plot{k}(:)'
        if strcmp(scaling,'lin'); plot(PAR.time./1000,      VAR.Cl(:,i) ,ls{:},lw{:},'Color',color(i,:)); end;
        if strcmp(scaling,'log'); plot(PAR.time./1000,log10(VAR.Cl(:,i)),ls{:},lw{:},'Color',color(i,:)); end;
    end
    % add equilibrium liquid composition if specified
    if PLOT.EQ_plot
        for i = PLOT.R_plot{k}(:)'
            if strcmp(scaling,'lin'); plot(PAR.time./1000,      VAR.ClEq(:,i) ,':',lw{1},0.5,'Color',color(i,:)); end;
            if strcmp(scaling,'log'); plot(PAR.time./1000,log10(VAR.ClEq(:,i)),':',lw{1},0.5,'Color',color(i,:)); end;
        end
    end
    % add fractional liquid composition if specified
    if PLOT.FR_plot
        for i = PLOT.R_plot{k}(:)'
            if strcmp(scaling,'lin'); plot(PAR.time./1000,      VAR.ClFr(:,i) ,'-.',lw{1},0.5,'Color',color(i,:)); end;
            if strcmp(scaling,'log'); plot(PAR.time./1000,log10(VAR.ClFr(:,i)),'-.',lw{1},0.5,'Color',color(i,:)); end;
        end
    end
    box on; grid on;
    set(gca,'xlim',[0,PAR.tend./1000],fs{:},tl{:});
    leg = legend(PAR.CompStr{PLOT.C_plot{k}(:)'}); set(leg,loc{:},ort{:},fs{:},tx{:});
    xlabel('Time [kyr]',fs{:},tx{:});
    ylabel('Melt composition [wt \%]',fs{:},tx{:});
    drawnow;
    if ~holdfig; hold off; end
end

%***  plot bulk composition with time
for k = 1:size(PLOT.CB_plot,2)
    figure(f); f=f+1;
    if ~holdfig; clf; end; hold on;
    for i = PLOT.CB_plot{k}(:)'
        if strcmp(scaling,'lin'); plot(PAR.time./1000,      VAR.C(:,i) ,ls{:},lw{:},'Color',color(i,:)); end;
        if strcmp(scaling,'log'); plot(PAR.time./1000,log10(VAR.C(:,i)),ls{:},lw{:},'Color',color(i,:)); end;
    end
    box on; grid on;
    set(gca,'xlim',[0,PAR.tend./1000],fs{:},tl{:});
    leg = legend(PAR.CompStr{PLOT.CB_plot{k}(:)'}); set(leg,loc{:},ort{:},fs{:},tx{:});
    xlabel('Time [kyr]',fs{:},tx{:});
    ylabel('Bulk composition [wt \%]',fs{:},tx{:});
    drawnow;
    if ~holdfig; hold off; end
end

end


%*****  define rates of change for equilibration by reaction  *************

function  [dYdt] = RatesOfChange(~,Y,PAR)

%***  read solution from ODE solver into data structure
VAR.P   =  Y(1);
VAR.T   =  Y(2).*1000;  % scale magnitude relative to other variables
VAR.f   =  Y(3);
VAR.Cs  =  Y(3+       (1:PAR.nc))';
VAR.Cl  =  Y(3+PAR.nc+(1:PAR.nc))';
VAR.C   =  repmat(VAR.f,1,PAR.nc).*VAR.Cl + (1-repmat(VAR.f,1,PAR.nc)).*VAR.Cs;
dYdt    =  zeros(size(Y));

%***  call R_DMC for new equilibrium f,C_s^i,C_l^i
[~  ,PAR]  =  R_DMC_Equilibrium(VAR,PAR,'K');
[VAReq,~]  =  R_DMC_Equilibrium(VAR,PAR,'E');

%***  calculate reaction rates Gamma^i
Gamma     =  R_DMC_ReactionRates(VAR,VAReq,PAR);
GammaNet  =  sum(Gamma);

%***  safeguard: fRate=0 at lower/upper bounds of melt fraction
PAR.fRate(VAR.f<=1e-6||(1-VAR.f)<=1e-6)  =  0; 

%***  set rate of change for pressure [GPa / yr]
dYdt(1)  =  PAR.PRate;

%***  set rate of change for temperature T [deg C / yr]
dYdt(2)  =  PAR.TRate + (PAR.alpha.*VAR.T.*PAR.PRate.*1e9 - sum(Gamma.*PAR.L))./PAR.rho./PAR.cp;
dYdt(2)  =  dYdt(2)./1000;  % scale magnitude relative to other variables

%***  set rate of change for melt fraction f [wt / yr]
dYdt(3)  =  PAR.fRate + GammaNet/PAR.rho;

%***  set rate of change for solid composition [wt % / yr]
for i = 1:PAR.nc;
    dYdt(3+i)         =  -(Gamma(i) - VAR.Cs(i).*GammaNet)./(max(1e-6,1-VAR.f).*PAR.rho);
    dYdt(3+PAR.nc+i)  =   (Gamma(i) - VAR.Cl(i).*GammaNet)./(max(1e-6,  VAR.f).*PAR.rho);
end

%***  update current reaction rates to pass as output
dYdt(end-PAR.nc+1:end)  =  (Gamma'-Y(end-PAR.nc+1:end))./(PAR.tend./PAR.nt);

end
