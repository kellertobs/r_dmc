%**************************************************************************
%*****  R_DMC PHASE DIAGRAM PLOTTING ROUTINE  *****************************
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
% FILENAME  R_DMC_PhaseDiagrams.m
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
% HELP:  [VAR,PAR] = R_DMC_PhaseDiagrams(PAR,PLOT)
% INPUT/OUTPUT: 
% - Input parameters structure PAR:
%   T0    melting point at P=0  [deg C  ]  1-by-nc row vector
%   A     linear P-coeff        [K/GPa  ]  "
%   B     quadratic P-coeff     [K/GPa^2]  "
%   L     latent heat           [J/kg   ]  " 
%   r     tuning parameter      [J/K/kg ]  "
%   nc    number of components  [1      ]  integer
%   np    number of gridpoints  [1      ]  integer
%   Pref  reference pressure    [GPa    ]  scalar
%   Pmin  minimum pressure      [GPa    ]  scalar
%   Pmax  maximum pressure      [GPa    ]  scalar
%   Tref  reference temperature [deg C  ]  scalar
%   Tmin  minimum temperature   [deg C  ]  scalar
%   Tmax  maximum temperature   [deg C  ]  scalar
%   Cref  reference bulk comp.  [wt %   ]  1-by-nc row vector
%   alpha thermal expansivity   [1/K    ]  scalar
%   cp    specific heat         [J/kg/K ]  scalar
%   rho   density               [kg/m3  ]  scalar
%   g     gravity               [m/s2   ]  scalar
% - Input plotting control structure PLOT:
%   hold     hold previous figures   [true/false ]  logical
%   LS       plotting line style     ['-',':',etc]  string
%   scale    choice of plot scaling  ['lin'/'log']  string
%   SL_plot  plot solidus/liquidus   [1/0]          integer
%   FT_plot  plot f vs T             [1/0]          integer
%   FP_plot  plot f vs P             [1/0]          integer
%   TP_plot  plot T_m^i vs P    e.g. {[1,2],[3,4]}  cell array
%   KP_plot  plot K^i vs P      e.g. {[1,2],[3,4]}  cell array
%   KT_plot  plot K^i vs T      e.g. {[1,2],[3,4]}  cell array
%   CT_plot  plot C^i vs T      e.g. {[1,2],[3,4]}  cell array
%   CP_plot  plot C^i vs P      e.g. {[1,2],[3,4]}  cell array
%   TC_bin_plot  plot binary phase loops   e.g. {[1,2],[3,4]}      cell array
%   TC_ter_plot  plot ternary phase loops  e.g. {[1,2,3],[2,3,4]}  cell array
% - Output loaded into VAR, PAR structures.


function  [VAR,PAR] = R_DMC_PhaseDiagrams(PAR,PLOT)

%*****  set up phase diagram plotting for specified input  ****************

%*** normalize input C_ref to 1
PAR.Cref  =  PAR.Cref./repmat(sum(PAR.Cref),1,PAR.nc);

%***  set up P_GPa vector for plotting
PAR.P_GPa  =  linspace(PAR.Pmin,PAR.Pmax,PAR.np)';

%***  set up P_GPa vector for plotting
PAR.T_degC  =  linspace(PAR.Tmin,PAR.Tmax,PAR.np)';

%***  set up T_adiabat vector for plotting
PAR.T_adiabat  =  PAR.Tref.*exp(PAR.P_GPa.*1e9.*PAR.alpha./PAR.rho./PAR.cp);

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
color    =  [flipud(R_DMC_Colormap(4:12:end,:));flipud(R_DMC_Colormap(8:12:end,:))];
set(0,'DefaultFigureColor',[1,1,1]);
set(0,'DefaultAxesLineWidth',1.5);
set(0,'DefaultFigureRenderer','zbuffer');

%***  plot component melting points: T_m^i vs P
if (~isempty(PLOT.TP_plot))
    % prepare solution variable structure VAR
    VAR.P  =  PAR.P_GPa;
    VAR.T  =  ones(size(VAR.P)).*PAR.Tref;
    VAR.C  =  repmat(PAR.Cref,size(VAR.P));
    % get component melting points T_m^i
    [~,PAR]  =  R_DMC_Equilibrium(VAR,PAR,'K');
    % create specified plots
    for k = 1:size(PLOT.TP_plot,2)
        figure(f); f=f+1;
        if ~holdfig; clf; end; hold on;
        for i = PLOT.TP_plot{k}(:)'
            plot(PAR.Tm(:,i),VAR.P,ls{:},lw{:},'Color',color(i,:));
        end
        box on; axis ij tight; grid on;
        set(gca,fs{:},tl{:});
        leg = legend(PAR.CompStr{PLOT.TP_plot{k}(:)'}); set(leg,loc{:},ort{:},fs{:},tx{:});
        xlabel('Melting Point [deg C]',fs{:},tx{:});
        ylabel('Pressure [GPa]',fs{:},tx{:});
        drawnow;
        if ~holdfig; hold off; end
    end
end

%***  plot component partition coefficients K^i vs T at P
if (~isempty(PLOT.KT_plot)) 
    % prepare solution variable structure VAR
    VAR.T  =  PAR.T_degC;
    VAR.P  =  ones(size(VAR.T)).*PAR.Pref;
    VAR.C  =  repmat(PAR.Cref,size(VAR.P));
    % get component partition coefficients K^i
    [~,PAR]  =  R_DMC_Equilibrium(VAR,PAR,'K');
    % create specified plots
    for k = 1:size(PLOT.KT_plot,2)
        figure(f); f=f+1;
        if ~holdfig; clf; end; hold on;
        for i = PLOT.KT_plot{k}(:)'
            plot(VAR.T,log10(PAR.K(:,i)),ls{:},lw{:},'Color',color(i,:));
        end
        box on; axis xy tight; grid on;
        set(gca,fs{:},tl{:});
        leg = legend(PAR.CompStr{PLOT.KT_plot{k}(:)'}); set(leg,loc{:},ort{:},fs{:},tx{:});
        xlabel('Temperature [deg C]',fs{:},tx{:});
        ylabel('Log Partition Coeff',fs{:},tx{:});
        drawnow;
        if ~holdfig; hold off; end
    end
end

%***  plot K^i vs P at T_adiabat
if (~isempty(PLOT.KP_plot))
    % prepare solution variable structure VAR
    VAR.P  =  PAR.P_GPa;
    VAR.T  =  PAR.T_adiabat;
    VAR.C  =  repmat(PAR.Cref,size(VAR.P));
    % get component partition coefficients K^i
    [~,PAR]  =  R_DMC_Equilibrium(VAR,PAR,'K');
    % create specified plots
    for k = 1:size(PLOT.KP_plot,2)
        figure(f); f=f+1;
        if ~holdfig; clf; end; hold on;
        for i = PLOT.KP_plot{k}(:)'
            plot(log10(PAR.K(:,i)),VAR.P,ls{:},lw{:},'Color',color(i,:));
        end
        box on; axis ij tight; grid on;
        set(gca,fs{:},tl{:});
        leg = legend(PAR.CompStr{PLOT.KP_plot{k}(:)'}); set(leg,loc{:},ort{:},fs{:},tx{:});
        xlabel('Log Partition Coeff',fs{:},tx{:});
        ylabel('Pressure [GPa]',fs{:},tx{:});
        drawnow;
        if ~holdfig; hold off; end
    end
end

%***  plot T_sol, T_liq vs P
if (~isempty(PLOT.SL_plot))
    % prepare solution variable structure VAR
    VAR.P  =  PAR.P_GPa;
    VAR.T  =  zeros(PAR.np,1);
    VAR.C  =  repmat(PAR.Cref,size(VAR.P));
    % get solidus and liquidus curves Tsol, Tliq
    [~,PAR]  =  R_DMC_Equilibrium(VAR,PAR,'T');
    % create specified plot
    figure(f); f=f+1;
    if ~holdfig; clf; end;  hold on;
    plot(PAR.Tsol,VAR.P,ls{:},lw{:},'Color',color(5,:));
    plot(PAR.Tliq,VAR.P,ls{:},lw{:},'Color',color(1,:));
    plot(PAR.T_adiabat,VAR.P,'-.',lw{:},'Color','k');
    box on; axis ij tight; grid on;
    set(gca,fs{:},tl{:});
    leg = legend('Solidus','Liquidus','Adiabat'); set(leg,loc{:},ort{:},fs{:},tx{:});
    xlabel('Temperature [deg C]',fs{:},tx{:});
    ylabel('Pressure [GPa]',fs{:},tx{:});
    drawnow;
    if ~holdfig; hold off; end
end

%***  plot binary TC phase diagram at Pref
if (~isempty(PLOT.TC_bin_plot))
    for k = 1:size(PLOT.TC_bin_plot,2)
        % get component pair for binary phase loop
        c1  =  PLOT.TC_bin_plot{k}(1);
        c2  =  PLOT.TC_bin_plot{k}(2);
        % prepare solution variable structure VAR
        VAR.C        =  zeros(PAR.np,PAR.nc);
        VAR.C(:,c2)  =  linspace(0,1,PAR.np);
        VAR.C(:,c1)  =  1-VAR.C(:,c2);
        VAR.P        =  ones(size(VAR.C(:,1))).*PAR.Pref;
        VAR.T        =  ones(size(VAR.C(:,1))).*PAR.Tref;
        % get solidus and liquidus curves Tsol, Tliq
        [~,PAR]  =  R_DMC_Equilibrium(VAR,PAR,'T');
        % create specified plots
        figure(f); f=f+1;
        if ~holdfig; clf; end; hold on;
        i = PLOT.TC_bin_plot{k}(2);
        plot(VAR.C(:,i).*100,PAR.Tsol,ls{:},lw{:},'Color',color(i,:));
        plot(VAR.C(:,i).*100,PAR.Tliq,ls{:},lw{:},'Color',color(i,:));
        box on; axis xy tight; grid on;
        set(gca,fs{:},tl{:});
        xlabel([PAR.CompStr{c1},' \& ',PAR.CompStr{c2},' [wt \%]'],fs{:},tx{:});
        ylabel('Temperature [deg C]',fs{:},tx{:});
        drawnow;
        if ~holdfig; hold off; end
    end
end

%***  plot ternary TC phase diagram at Pref
if (~isempty(PLOT.TC_ter_plot))    
    for k = 1:size(PLOT.TC_ter_plot,2)
        % get three components for ternary phase loop
        c1  =  PLOT.TC_ter_plot{k}(1);
        c2  =  PLOT.TC_ter_plot{k}(2);
        c3  =  PLOT.TC_ter_plot{k}(3);
        % prepare solution variable structure VAR
        C3  =  linspace(0,1,PAR.np);
        C2  =  linspace(0,1,PAR.np);
        [C2,C3]  =  meshgrid(C2,C3);
        C2  =  C2(:);  C3  =  C3(:);  C1  =  1-C2-C3;
        VAR.C  =  zeros(PAR.np*PAR.np,PAR.nc);
        VAR.C(:,c1) = C1;  VAR.C(:,c2) = C2;  VAR.C(:,c3) = C3;
        VAR.P  =  ones(size(VAR.C(:,1))).*PAR.Pref;
        VAR.T  =  ones(size(VAR.C(:,1))).*PAR.Tref;
        % get solidus and liquidus curves Tsol, Tliq
        [~,PAR]  =  R_DMC_Equilibrium(VAR,PAR,'T');
        % exclude points of invalid composition
        PAR.Tsol(VAR.C(:,c2)+VAR.C(:,c3)>1)  =  nan;
        PAR.Tliq(VAR.C(:,c2)+VAR.C(:,c3)>1)  =  nan;
        % plot ternary axes
        figure(f); f=f+1; clf
        hold on;
        set(gca,'visible','off');
        patch('xdata',[0,1,0.5,0],'ydata',[0,0,sin(1/3*pi),0],'edgecolor',[1,1,1],'facecolor',[1,1,1],'handlevisibility','off');
        plot([0,1,0.5,0],[0,0,sin(1/3*pi),0],'color',[0,0,0],'handlevisibility','off');
        % plot grid on ternary axes
        m = 5;
        grids = linspace(0,100,m+1);
        grids = grids(1:end-1);
        labels = num2str(grids(2:end)');
        [x3, y3] = terncoords(100-grids, grids, zeros(size(grids)));
        [x2, y2] = terncoords(grids, zeros(size(grids)), 100-grids);
        [x1, y1] = terncoords(zeros(size(grids)), 100-grids, grids);
        n = m-1;
        for i = 1:n
            plot([x1(i+1),x2(n-i+2)],[y1(i+1),y2(n-i+2)],'k:','handlevisibility','off');
            plot([x2(i+1),x3(n-i+2)],[y2(i+1),y3(n-i+2)],'k:','handlevisibility','off');
            plot([x3(i+1),x1(n-i+2)],[y3(i+1),y1(n-i+2)],'k:','handlevisibility','off');
        end
        % plot T_sol, T_liq on ternary axes
        [X,Y]  =  terncoords(VAR.C(:,c2),VAR.C(:,c3),VAR.C(:,c1));
        tri    =  simpletri(PAR.np);
        trisurf(tri,X,Y,PAR.Tsol);
        trisurf(tri,X,Y,PAR.Tliq);
        shading interp;
        view([60,35]);
        pos = get(gca,'position'); pos(1) = pos(1)+0.09; pos(2) = pos(2)+.18; set(gca,'position',pos);
        cb = colorbar; set(cb,fs{:},tl{:},'Position',[0.1,0.16,0.8,0.05],ort{:}); colormap(R_DMC_Colormap);
        text(0.054,-0.3,'Solidus \& Liquidus T [deg C]','Units','normalized',fs{:},tx{:});
        text(x3(2:end)+0.09,y3(2:end)-0.04,labels,fs{:},tx{:});
        text(x2(2:end)+0.04,y2(2:end)-0.12,labels,fs{:},tx{:});
        text(x1(2:end)-0.09,y1(2:end)+0.02,labels,fs{:},tx{:});
        text( 0.3686, 0.8845,PAR.CompStr{c3},fs{:},tx{:});
        text( 1.0302,-0.1375,PAR.CompStr{c2},fs{:},tx{:});
        text(-0.1929,-0.0840,PAR.CompStr{c1},fs{:},tx{:});
        drawnow;
    end
end

%***  plot evolution of f,Cs,Cl with T
if (~isempty(PLOT.FT_plot) || ~isempty(PLOT.CT_plot))
    % prepare solution variable structure VAR
    VAR.T  =  PAR.T_degC;
    VAR.P  =  ones(size(VAR.T)).*PAR.Pref;
    VAR.f  =  zeros(size(VAR.T));
    VAR.C  =  repmat(PAR.Cref,size(VAR.P));
    % get equilibrium melt and component fractions
    [VAR,~]  =  R_DMC_Equilibrium(VAR,PAR,'E');
    % exclude solid/liquid compositions where there is no solid/liquid phase
    VAR.Cs(1-VAR.f<1e-6,:)  =  nan;  PAR.CsFract(1-VAR.f<1e-6,:)  =  nan;
    VAR.Cl(  VAR.f<1e-6,:)  =  nan;  PAR.ClFract(  VAR.f<1e-6,:)  =  nan;
    % plot evolution of f with T
    if (~isempty(PLOT.FT_plot))
        figure(f); f=f+1;
        if ~holdfig; clf; end; hold on;
        if strcmp(scaling,'lin'); plot(VAR.T,               VAR.f .*100 ,ls{:},lw{:},'Color',color(1,:)); end;
        if strcmp(scaling,'log'); plot(VAR.T,log10(max(1e-6,VAR.f).*100),ls{:},lw{:},'Color',color(1,:)); end;
        box on; axis xy tight; grid on;
        set(gca,fs{:},tl{:});
        xlabel('Temperature [deg C]',fs{:},tx{:});
        ylabel('Melt Fraction [wt \%]',fs{:},tx{:});
        drawnow;
        if ~holdfig; hold off; end
    end
    % plot evolution of C_s^i with T
    for k = 1:size(PLOT.CT_plot,2)
        figure(f); f=f+1;
        if ~holdfig; clf; end; hold on;
        for i = PLOT.CT_plot{k}(:)'
            if strcmp(scaling,'lin'); plot(VAR.T,      VAR.Cs(:,i).*100 ,ls{:},lw{:},'Color',color(i,:)); end
            if strcmp(scaling,'log'); plot(VAR.T,log10(VAR.Cs(:,i).*100),ls{:},lw{:},'Color',color(i,:)); end
        end
        box on; axis xy tight; grid on;
        set(gca,fs{:},tl{:});
        leg = legend(PAR.CompStr{PLOT.CT_plot{k}(:)'}); set(leg,loc{:},ort{:},fs{:},tx{:});
        xlabel('Temperature [deg C]',fs{:},tx{:});
        ylabel('Rock Composition [wt \%]',fs{:},tx{:});
        drawnow;
        if ~holdfig; hold off; end
    end
    % plot evolution of C_l^i with T
    if (~isempty(PLOT.CT_plot))
        for k = 1:size(PLOT.CT_plot,2)
            figure(f); f=f+1;
            if ~holdfig; clf; end; hold on;
            for i = PLOT.CT_plot{k}(:)'
                if strcmp(scaling,'lin'); plot(VAR.T,      VAR.Cl(:,i).*100 ,ls{:},lw{:},'Color',color(i,:)); end
                if strcmp(scaling,'log'); plot(VAR.T,log10(VAR.Cl(:,i).*100),ls{:},lw{:},'Color',color(i,:)); end
            end
            box on; axis xy tight; grid on;
            set(gca,fs{:},tl{:});
            leg = legend(PAR.CompStr{PLOT.CT_plot{k}(:)'}); set(leg,loc{:},ort{:},fs{:},tx{:});
            xlabel('Temperature [deg C]',fs{:},tx{:});
            ylabel('Melt Composition [wt \%]',fs{:},tx{:});
            drawnow;
            if ~holdfig; hold off; end
        end
    end
end

%***  plot evolution of f,Cs,Cl with P
if (~isempty(PLOT.FP_plot) || ~isempty(PLOT.CP_plot))
    % prepare solution variable structure VAR
    VAR.P  =  PAR.P_GPa;
    VAR.T  =  PAR.T_adiabat;
    VAR.f  =  zeros(size(VAR.P));
    VAR.C  =  repmat(PAR.Cref,size(VAR.P));
    % get equilibrium melt and component fractions
    [VAR,~]  =  R_DMC_Equilibrium(VAR,PAR,'E');
    % exclude solid/liquid compositions where there is no solid/liquid phase
    VAR.Cs(1-VAR.f<1e-6,:)  =  nan;  PAR.CsFract(1-VAR.f<1e-6,:)  =  nan;
    VAR.Cl(  VAR.f<1e-6,:)  =  nan;  PAR.ClFract(  VAR.f<1e-6,:)  =  nan;
    % plot evolution of f with P
    if (~isempty(PLOT.FP_plot))
        figure(f); f=f+1;
        if ~holdfig; clf; end; hold on;
        if strcmp(scaling,'lin'); plot(               VAR.f .*100 ,VAR.P,ls{:},lw{:},'Color',color(1,:)); end;
        if strcmp(scaling,'log'); plot(log10(max(1e-6,VAR.f).*100),VAR.P,ls{:},lw{:},'Color',color(1,:)); end;
        box on; axis ij tight; grid on;
        set(gca,fs{:},tl{:});
        xlabel('Melt Fraction [wt \%]',fs{:},tx{:});
        ylabel('Pressure [GPa]',fs{:},tx{:});
        drawnow;
        if ~holdfig; hold off; end
    end
    % plot evolution of C_l^i with P
    for k = 1:size(PLOT.CP_plot,2)
        figure(f); f=f+1;
        if ~holdfig; clf; end; hold on;
        for i = PLOT.CP_plot{k}(:)'
            if strcmp(scaling,'lin'); plot(      VAR.Cs(:,i).*100 ,VAR.P,ls{:},lw{:},'Color',color(i,:)); end;
            if strcmp(scaling,'log'); plot(log10(VAR.Cs(:,i).*100),VAR.P,ls{:},lw{:},'Color',color(i,:)); end;
        end
        box on; axis ij tight; grid on;
        set(gca,fs{:},tl{:});
        leg = legend(PAR.CompStr{PLOT.CP_plot{k}(:)'}); set(leg,loc{:},ort{:},fs{:},tx{:});
        xlabel('Rock Composition [wt \%]',fs{:},tx{:});
        ylabel('Pressure [GPa]',fs{:},tx{:});
        drawnow;
        if ~holdfig; hold off; end
    end
    % plot evolution of C_l^i with P
    if (~isempty(PLOT.CP_plot))
        for k = 1:size(PLOT.CP_plot,2)
            figure(f); f=f+1;
            if ~holdfig; clf; end; hold on;
            for i = PLOT.CP_plot{k}(:)'
                if strcmp(scaling,'lin'); plot(      VAR.Cl(:,i).*100 ,VAR.P,ls{:},lw{:},'Color',color(i,:)); end;
                if strcmp(scaling,'log'); plot(log10(VAR.Cl(:,i).*100),VAR.P,ls{:},lw{:},'Color',color(i,:)); end;
            end
            box on; axis ij tight; grid on;
            set(gca,fs{:},tl{:});
            leg = legend(PAR.CompStr{PLOT.CP_plot{k}(:)'}); set(leg,loc{:},ort{:},fs{:},tx{:});
            xlabel('Melt Composition [wt \%]',fs{:},tx{:});
            ylabel('Pressure [GPa]',fs{:},tx{:});
            drawnow;
            if ~holdfig; hold off; end
        end
    end 
end

end


% The following are utility functions for plotting of ternary phase diagrams, 
% adapted from TERNPLOT package by Carl Sandrock, 2012:
% www.mathworks.com/matlabcentral/fileexchange/2299-alchemyst-ternplot

function tri = simpletri(N)
[row, col] = meshgrid(1:N-1);

bl = sub2ind([N, N], row, col);
bl = bl(:);
br = bl + 1;
tl = bl + N;
tr = tl + 1;

tri = [tl bl br; tl tr br];
end

function [x,y] = terncoords(c1,c2,c3)
csum  =  c1 + c2 + c3;
c1  =  c1./csum;  c2 = c2./csum;  c3 = c3./csum;

y = c2     .* sin(deg2rad(60));
x = c1 + y .* cot(deg2rad(60));
end

function rad = deg2rad(deg)
rad = deg / 180 * pi;
end