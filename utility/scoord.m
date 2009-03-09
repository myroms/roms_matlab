function [z,s,C]=scoord(h, x, y, Vtransform, Vstretching, ...
                        theta_s, theta_b, hc,             ...
			N, kgrid, column, index, plt);
%
% SCOORD:  Compute and plot ROMS vertical stretched coordinates
%
% [z,s,C]=scoord(h, x, y, Vtransform, Vstretching, theta_s, theta_b, ...
%                hc, N, kgrid, column, index, plt)
%
% Given a batymetry (h) and terrain-following stretching parameters,
% this function computes the depths of RHO- or W-points for a vertical
% grid section along columns (ETA-axis) or rows (XI-axis). Check the
% following link for details:
%
%    https://www.myroms.org/wiki/index.php/Vertical_S-coordinate
%
% On Input:
%
%    h             Bottom depth 2D array, h(1:Lp,1:Mp), m, positive
%    x             X-coordinate 2D array, x(1:Lp,1:Mp), m or degrees_east
%    y             Y-coordinate 2D array, y(1:Lp,1:Mp), m or degrees_north
%    Vtransform    Vertical transformation equation:
%                    Vtransform = 1,   original transformation
%
%                      z(x,y,s,t)=Zo(x,y,s)+zeta(x,y,t)*[1+Zo(x,y,s)/h(x,y)]
%
%                      Zo(x,y,s)=hc*s+[h(x,y)-hc]*C(s)
%
%                    Vtransform = 2,   new transformation
%
%                      z(x,y,s,t)=zeta(x,y,t)+[zeta(x,y,t)+h(x,y)]*Zo(x,y,s)
%
%                       Zo(x,y,s)=[hc*s(k)+h(x,y)*C(k)]/[hc+h(x,y)]
%    Vstretching   Vertical stretching function:
%                    Vstretching = 1,  original (Song and Haidvogel, 1994)
%                    Vstretching = 2,  A. Shchepetkin (UCLA-ROMS)
%                    Vstretching = 3,  R. Geyer BBL refinement
%    theta_s       S-coordinate surface control parameter (scalar)
%    theta_b       S-coordinate bottom control parameter (scalar)
%    hc            Width (m) of surface or bottom boundary layer in which
%                    higher vertical resolution is required during
%                    stretching (scalar)
%    N             Number of vertical levels (scalar)
%    kgrid         Depth grid type logical switch:
%                    kgrid = 0,        depths of RHO-points
%                    kgrid = 1,        depths of W-points
%    column        Grid direction logical switch:
%                    column = 1,       column section
%                    column = 0,       row section
%    index         Column or row to compute (scalar)
%                    if column = 1,    then   1 <= index <= Lp
%                    if column = 0,    then   1 <= index <= Mp
%    plt           Switch to plot scoordinate (scalar):
%                    plt = 0,          do not plot
%                    plt = 1,          plot
%
% On Output:
%
%    z             Depths (m) of RHO- or W-points (matrix)
%    s             S-coordinate independent variable, [-1 < s < 0] at
%                    vertical RHO- or W-points (vector)
%    C             Set of S-curves used to stretch the vertical coordinate
%                    lines that follow the topography at vertical RHO- or
%                    W-points (vector)
%

% svn $Id$ 
%===========================================================================%
%  Copyright (c) 2002-2009 The ROMS/TOMS Group                              %
%    Licensed under a MIT/X style license                                   %
%    See License_ROMS.txt                           Hernan G. Arango        %
%===========================================================================%


%----------------------------------------------------------------------------
%  Set several parameters.
%----------------------------------------------------------------------------

if (hc > min(min(h))),
  disp([setstr(7),'*** Error:  SCOORD - critical depth exceeds minimum' ...
	' bathymetry value.',setstr(7)]);
  disp([setstr(7),'                     hc   = ',num2str(hc),setstr(7)]);
  disp([setstr(7),'                     hmax = ',num2str(min(min(h))), ...
	setstr(7)]);
  return
end,

if (Vtransform < 1 | Vtransform > 2),
  disp([setstr(7),'*** Error:  SCOORD - Illegal parameter Vtransform = ' ...
	num2str(Vtransfrom), setstr(7)]);
end,

if (Vstretching < 1 | Vstretching > 3),
  disp([setstr(7),'*** Error:  SCOORD - Illegal parameter Vstretching = ' ...
	num2str(Vstretching), setstr(7)]);
end,

Np=N+1;
[Lp Mp]=size(h);
hmin=min(min(h));
hmax=max(max(h));

%----------------------------------------------------------------------------
% Test input to see if it's in an acceptable form.
%----------------------------------------------------------------------------

if (nargin < 12),
  disp(' ');
  disp([setstr(7),'*** Error:  SCOORD - too few arguments.',setstr(7)]);
  disp([setstr(7),'                     number of supplied arguments: ',...
       num2str(nargin),setstr(7)]);
  disp([setstr(7),'                     number of required arguments: 8',...
       setstr(7)]);
  disp(' ');
  return
end,

if (column),

  if (index < 1 | index > Lp),
    disp(' ');
    disp([setstr(7),'*** Error:  SCOORD - illegal column index.',setstr(7)]);
    disp([setstr(7),'                     valid range:  1 <= index <= ',...
         num2str(Lp),setstr(7)]);
    disp(' ');
    return
  end,

else,

  if (index < 1 | index > Mp),
    disp(' ');
    disp([setstr(7),'*** Error:  SCOORD - illegal row index.',setstr(7)]);
    disp([setstr(7),'                     valid range:  1 <= index <= ',...
         num2str(Mp),setstr(7)]);
    disp(' ');
    return
  end,

end,

%----------------------------------------------------------------------------
% Define S-Curves at vertical RHO- or W-points (-1 < s < 0).
%----------------------------------------------------------------------------

% Original vertical stretching function (Song and Haidvogel, 1994).

if (Vstretching == 1),

  ds=1.0/N;
  if (kgrid == 1),
    Nlev=Np;
    lev=0:N;
    s=(lev-N).*ds;
  else,
    Nlev=N;
    lev=[1:N]-0.5;
    s=(lev-N).*ds;
  end,
  if (theta_s > 0),
    Ptheta=sinh(theta_s.*s)./sinh(theta_s);
    Rtheta=tanh(theta_s.*(s+0.5))./(2.0*tanh(0.5*theta_s))-0.5;
    C=(1.0-theta_b).*Ptheta+theta_b.*Rtheta;
  else,
    C=s;
  end,

% A. Shchepetkin (UCLA-ROMS) vertical stretching function.

elseif (Vstretching == 2),

  alfa=1.0;
  beta=1.0;
  ds=1.0/N;
  if (kgrid == 1),
    Nlev=Np;
    lev=0:N;
    s=(lev-N).*ds;
  else
    Nlev=N;
    lev=[1:N]-0.5;
    s=(lev-N).*ds;
  end,
  if (theta_s > 0),
    Csur=(1.0-cosh(theta_s.*s))/(cosh(theta_s)-1.0);
    if (theta_b > 0),
      Cbot=-1.0+sinh(theta_b*(s+1.0))/sinh(theta_b);
      weigth=(s+1.0).^alfa.*(1.0+(alfa/beta).*(1.0-(s+1.0).^beta));
      C=weigth.*Csur+(1.0-weigth).*Cbot;
    else
      C=Csur;
    end,
  else,
    C=s;
  end,
  
%  R. Geyer BBL vertical stretching function.

elseif (Vstretching == 3),

  ds=1.0/N;
  if (kgrid == 1),
    Nlev=Np;
    lev=0:N;
    s=(lev-N).*ds;
  else,
    Nlev=N;
    lev=[1:N]-0.5;
    s=(lev-N).*ds;
  end,
  if (theta_s > 0),
     exp_s=theta_s;      %  surface stretching exponent
     exp_b=theta_b;      %  bottom  stretching exponent
     alpha=3;            %  scale factor for all hyperbolic functions
    Cbot=log(cosh(alpha*(s+1).^exp_b))/log(cosh(alpha))-1;
    Csur=-log(cosh(alpha*abs(s).^exp_s))/log(cosh(alpha));
    weight=(1-tanh( alpha*(s+.5)))/2;
    C=weight.*Cbot+(1-weight).*Csur;
  else,
    C=s;
  end,

end,

% Report S-coordinate parameters.

report = 1;
if report,
  disp(' ');
  if (Vstretching == 1),
    disp(['Vstretching = ',num2str(Vstretching), '   Song and Haidvogel (1994)']);
  elseif (Vstretching == 2),
    disp(['Vstretching = ',num2str(Vstretching), '   Shchepetkin (2005)']);
  elseif (Vstretching == 3),
    disp(['Vstretching = ',num2str(Vstretching), '   Geyer (2009), BBL']);
  end,
  if (Vtransform == 1),
    disp(['Vtransform  = ',num2str(Vtransform), '   original ROMS']);
  elseif (Vtransform == 2),
    disp(['Vtransform  = ',num2str(Vtransform), '   ROMS-UCLA']);
  end,
  if (kgrid == 1)
    disp(['   kgrid    = ',num2str(kgrid), '   at W-points']);
  else,
    disp(['   kgrid    = ',num2str(kgrid), '   at RHO-points']);
  end,
  disp(['   theta_s  = ',num2str(theta_s)]);
  disp(['   theta_b  = ',num2str(theta_b)]);
  disp(['   hc       = ',num2str(hc)]);
  disp(['   hmin     = ',num2str(hmin)]);
  disp(['   hmax     = ',num2str(hmax)]);
  disp(' ');
  disp(' S-coordinate curves: k, s(k), C(k)')
  disp(' ');
  if (kgrid == 1),
    for k=Nlev:-1:1,
      disp(['    ', ...
	    sprintf('%3g',k-1     ), '   ', ...
	    sprintf('%20.12e',s(k)), '   ', ...
	    sprintf('%20.12e',C(k))]);
    end,
  else
    for k=Nlev:-1:1,
      disp(['    ', ...
	    sprintf('%3g',k       ), '   ', ...
	    sprintf('%20.12e',s(k)), '   ', ...
	    sprintf('%20.12e',C(k))]);
    end,
  end,
  disp(' ');
end,

%============================================================================
% Compute depths at requested grid section.  Assume zero free-surface.
%============================================================================

zeta=zeros(size(h));

%----------------------------------------------------------------------------
% Column section: section along ETA-axis.
%----------------------------------------------------------------------------

if (column),

  if (Vtransform == 1),

    z=zeros(Mp,Nlev);
    for k=1:Nlev,
      z0=hc.*(s(k)-C(k))+h(index,:)*C(k);
      z(:,k)=z0+zeta(index,:).*(1.0+z0/h(index,:));
    end,
   
  elseif (Vtransform == 2),
  
    z=zeros(Mp,Nlev);
    for k=1:Nlev,
      z0=(hc.*s(k)+C(k).*h(index,:))./(h(index,:)+hc);
      z(:,k)=zeta(index,:)+(zeta(index,:)+h(index,:)).*z0;
    end,

  end,

%----------------------------------------------------------------------------
% Row section: section along XI-axis.
%----------------------------------------------------------------------------

else,

  if (Vtransform == 1),

    z=zeros(Lp,Nlev);
    for k=1:Nlev,
      z0=hc.*(s(k)-C(k))+h(:,index)*C(k);
      z(:,k)=z0+zeta(:,index).*(1.0+z0./h(:,index));
    end,

  elseif (Vtransform == 2),
  
    z=zeros(Lp,Nlev);
    for k=1:Nlev,
      z0=(hc.*s(k)+C(k).*h(:,index))./(h(:,index)+hc);
      z(:,k)=zeta(:,index)+(zeta(:,index)+h(:,index)).*z0;
    end,
    
  end,
  
end,

%============================================================================
% Plot grid section.
%============================================================================

if nargin < 9
  plt = 1;
end

if (plt == 1),

  figure;
  
  if (column),

    eta=y(index,:);
    eta2=[eta(1) eta eta(Mp)];
    hs=-h(index,:);
    zmin=min(hs);
    hs=[zmin hs zmin];
    hold off;
    fill(eta2,hs,[0.6 0.7 0.6]);
    hold on;
    han=plot(eta',z);
    set(han,'color', [0.5 0.5 0.5]);

    set(gcf,'Units','Normalized',...
        'Position',[0.2 0.1 0.6 0.8],...
	'PaperUnits','Normalized',...
	'PaperPosition',[0.2 0.1 0.6 0.8]);
    set(gca,'xlim',[-Inf Inf],'ylim',[zmin 0]);
    if (kgrid == 0),
      title(['Grid (\rho-points) Section at  \xi = ',num2str(index)]);
    else
      title(['Grid (W-points) Section at  \xi = ',num2str(index)]);
    end,
    xlabel({['Vcoord = ', num2str(Vtransform),',',num2str(Vstretching), ...
	     '   \theta_s = ' num2str(theta_s), ...
             '   \theta_b = ' num2str(theta_b), ...
	     '   hc  = ' num2str(hc), ...
             '   N = ' num2str(N)],['y-axis']});
    ylabel('depth  (m)');

  else,

    xi=x(:,index)';
    xi2=[xi(1) xi xi(Lp)];
    hs=-h(:,index)';
    zmin=min(hs);
    hs=[zmin hs zmin];
    hold off;
    fill(xi2,hs,[0.6 0.7 0.6]);
    hold on;
    han=plot(xi,z);
    set(han,'color', [0.5 0.5 0.5]);

    set(gcf,'Units','Normalized',...
        'Position',[0.2 0.1 0.6 0.8],...
	'PaperUnits','Normalized',...
	'PaperPosition',[0.2 0.1 0.6 0.8]);
    set(gca,'xlim',[-Inf Inf],'ylim',[zmin 0]);
    if (kgrid == 0),
      title(['Grid (\rho-points) Section at  \eta = ',num2str(index)]);
    else
      title(['Grid (W-points) Section at  \eta = ',num2str(index)]);
    end
    xlabel({['Vcoord = ', num2str(Vtransform),',',num2str(Vstretching), ...
	     '   \theta_s = ' num2str(theta_s), ...
             '   \theta_b = ' num2str(theta_b), ...
	     '   hc  = ' num2str(hc), ...
             '   N = ' num2str(N)],['x-axis']});
    ylabel('depth  (m)');

  end,

end

return


