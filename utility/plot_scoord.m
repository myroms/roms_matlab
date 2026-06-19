function [z,s,C]=plot_scoord(G, kgrid, column, index, plt, Zzoom);
%
% plot_SCOORD:  Compute and plot ROMS vertical stretched coordinates
%
% [z,s,C]=plot_scoord(G, kgrid, column, index, plt, Zzoom)
%
% Given a ROMS full grid structure this function computes the
% depths of RHO- or W-points for a vertical grid section along columns
% (ETA-axis) or rows (XI-axis). Check the following link for details:
%
%    https://www.myroms.org/wiki/index.php/Vertical_S-coordinate
%
% This script is the same as "scoord.m" but with a full grid structure
% argument computed elsewhere:
%
%    G = get_roms_grid('roms_grd.nc', 'roms_his.nc', rec);
% or
%    G = get_roms_grid('roms_his.nc', 'roms_his.nc', rec);
%
% the record is optional to specify a non zero free-surface.
%
% On Input:
%
%    G             An existing ROMS full grid extructure (struc)
%                    It must contains vertical grid fields
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
%                    plt = 2,          plot 2 pannels with zoom
%    Zzoom         If plt=2, maximum depth of the zoom in upper pannel
%
% On Output:
%
%    z             Depths (m) of RHO- or W-points (matrix)
%

% svn $Id$
%===========================================================================%
%  Copyright (c) 2002-2026 The ROMS Group                                   %
%    Licensed under a MIT/X style license                                   %
%    See License_ROMS.md                            Hernan G. Arango        %
%===========================================================================%

z=[];

%------------------------------------------------------------------------
%  Set several parameters.
%------------------------------------------------------------------------

if (~isstruct(G))
  disp(' ');
  disp([setstr(7),'*** Error:  GET_SCOORD - not a ROMS structue.',    ...
        setstr(7)]);
end

[Lp Mp]=size(G.h);
hmin=min(min(G.h));
hmax=max(max(G.h));
havg=0.5*(hmax+hmin);

%------------------------------------------------------------------------
% Test input to see if it's in an acceptable form.
%------------------------------------------------------------------------

if (column)
  if (index < 1 | index > Lp)
    disp(' ');
    disp([setstr(7),'*** Error:  GET_SCOORD - illegal column index.', ...
          setstr(7)]);
    disp([setstr(7),'            valid range:  1 <= index <= ',       ...
          num2str(Lp),setstr(7)]);
    disp(' ');
    return
  end
else
  if (index < 1 | index > Mp)
    disp(' ');
    disp([setstr(7),'*** Error:  GET_SCOORD - illegal row index.',    ...
          setstr(7)]);
    disp([setstr(7),'            valid range:  1 <= index <= ',       ...
          num2str(Mp),setstr(7)]);
    disp(' ');
    return
  end
end

%----------------------------------------------------------------------------
% Compute vertical stretching function, C(k):
%----------------------------------------------------------------------------

[s,C]=stretching(G.Vstretching, G.theta_s, G.theta_b, G.hc, G.N, kgrid, 0);

if (kgrid == 1)
  Nlev=G.N+1;
else
  Nlev=G.N;
end

if (G.Vtransform == 1)

  for k=Nlev:-1:1,
    zhc(k)=G.hc*s(k);
    z1 (k)=zhc(k)+(hmin-G.hc)*C(k);
    z2 (k)=zhc(k)+(havg-G.hc)*C(k);
    z3 (k)=zhc(k)+(hmax-G.hc)*C(k);
  end,

elseif (G.Vtransform == 2)

  for k=Nlev:-1:1
    if (G.hc > hmax)
      zhc(k)=hmax*(G.hc*s(k)+hmax*C(k))/(G.hc+hmax);
    else
      zhc(k)=0.5*min(G.hc,hmax)*(s(k)+C(k));
    end
    z1 (k)=hmin*(G.hc*s(k)+hmin*C(k))/(G.hc+hmin);
    z2 (k)=havg*(G.hc*s(k)+havg*C(k))/(G.hc+havg);
    z3 (k)=hmax*(G.hc*s(k)+hmax*C(k))/(G.hc+hmax);
  end

end,

report=false;

if (report)
  disp(' ');
  if (G.Vtransform == 1)
    disp(['Vtransform  = ',num2str(G.Vtransform), '   original ROMS']);
  elseif (G.Vtransform == 2)
    disp(['Vtransform  = ',num2str(G.Vtransform), '   ROMS-UCLA']);
  end
  if (G.Vstretching == 1)
    disp(['Vstretching = ',num2str(G.Vstretching), '   Song and Haidvogel (1994)']);
  elseif (G.Vstretching == 2)
    disp(['Vstretching = ',num2str(G.Vstretching), '   Shchepetkin (2005)']);
  elseif (G.Vstretching == 3)
    disp(['Vstretching = ',num2str(G.Vstretching), '   Geyer (2009), BBL']);
  elseif (G.Vstretching == 4)
    disp(['Vstretching = ',num2str(G.Vstretching), '   Shchepetkin (2010)']);
  end,
  if (kgrid == 1)
    disp(['   kgrid    = ',num2str(kgrid), '   at vertical W-points']);
  else
    disp(['   kgrid    = ',num2str(kgrid), '   at vertical RHO-points']);
  end
  disp(['   theta_s  = ',num2str(G.theta_s)]);
  disp(['   theta_b  = ',num2str(G.theta_b)]);
  disp(['   hc       = ',num2str(G.hc)]);

  disp(' ');
  disp(' S-coordinate curves: ')
  disp(' ');
  disp([' level     S-coord    Cs-Curve   Z  at hmin      ', ...
        ' at hc    half way     at hmax']);
  disp(' ');

  if (kgrid == 1)
    for k=Nlev:-1:1
      disp(['   ',                         ...
            sprintf('%3i',k-1      ), ' ', ...
            sprintf('%12.7f',s(k)  ),      ...
            sprintf('%12.7f',C(k)  ),      ...
            sprintf('%12.3f',z1(k) ),      ...
            sprintf('%12.3f',zhc(k)),      ...
            sprintf('%12.3f',z2(k) ),      ...
            sprintf('%12.3f',z3(k) )]);
    end
  else
    for k=Nlev:-1:1
      disp(['   ',                         ...
            sprintf('%3i',k        ), ' ', ...
            sprintf('%12.7f',s(k)  ),      ...
            sprintf('%12.7f',C(k)  ),      ...
            sprintf('%12.3f',z1(k) ),      ...
            sprintf('%12.3f',zhc(k)),      ...
            sprintf('%12.3f',z2(k) ),      ...
            sprintf('%12.3f',z3(k) )]);
    end
  end
  disp(' ');

end

%============================================================================
% Compute depths at requested grid section.  Assume zero free-surface.
%============================================================================

zeta=zeros(size(G.h));

%----------------------------------------------------------------------------
% Column section: section along ETA-axis.
%----------------------------------------------------------------------------

if (column)

  if (G.Vtransform == 1)

    z=zeros(Mp,Nlev);
    for k=1:Nlev,
      z0=G.hc.*(s(k)-C(k))+G.h(index,:)*C(k);
      z(:,k)=z0+zeta(index,:).*(1.0+z0/G.h(index,:));
    end

  elseif (G.Vtransform == 2),

    z=zeros(Mp,Nlev);
    for k=1:Nlev
      z0=(G.hc.*s(k)+C(k).*G.h(index,:))./(G.h(index,:)+G.hc);
      z(:,k)=zeta(index,:)+(zeta(index,:)+G.h(index,:)).*z0;
    end

  end

%----------------------------------------------------------------------------
% Row section: section along XI-axis.
%----------------------------------------------------------------------------

else

  if (G.Vtransform == 1)

    z=zeros(Lp,Nlev);
    for k=1:Nlev
      z0=G.hc.*(s(k)-C(k))+G.h(:,index)*C(k);
      z(:,k)=z0+zeta(:,index).*(1.0+z0./G.h(:,index));
    end

  elseif (G.Vtransform == 2)

    z=zeros(Lp,Nlev);
    for k=1:Nlev,
      z0=(G.hc.*s(k)+C(k).*G.h(:,index))./(G.h(:,index)+G.hc);
      z(:,k)=zeta(:,index)+(zeta(:,index)+G.h(:,index)).*z0;
    end

  end

end

%========================================================================
% Plot grid section.
%========================================================================

if (nargin < 5)
  plt = 1;
end

if (plt > 0)

  figure;

  if (column)

    set(gcf,'Units','Normalized',           ...
        'Position',[0.2 0.1 0.6 0.8],       ...
        'PaperOrientation', 'landscape',    ...
        'PaperUnits','Normalized',          ...
        'PaperPosition',[0.2 0.1 0.6 0.8]);

    if (G.spherical)
      eta=G.lat_rho(index,:);
    else
      eta=G.y_rho(index,:);
    end
    eta2=[eta(1) eta eta(Mp)];
    hs=-G.h(index,:);
    zmin=min(hs);
    hs=[zmin hs zmin];

    if (plt == 2)
      h1=subplot(2,1,1);
      p1=get(h1,'pos');
      p1(2)=p1(2)+0.1;
      p1(4)=p1(4)-0.1;
      set (h1,'pos',p1);
    end

    hold off;
    fill(eta2,hs,[0.6 0.7 0.6]);
    hold on;
    han1=plot(eta',z);
    set(han1,'color', [0.5 0.5 0.5]);
    if (plt == 2)
      set(gca,'xlim',[-Inf Inf],'ylim',[-abs(Zzoom) 0]);
      ylabel('depth  (m)');
    else
      set(gca,'xlim',[-Inf Inf],'ylim',[zmin 0]);
    end

    if (kgrid == 0)
      title(['Grid (\rho-points) Section at  \xi = ',num2str(index)]);
    else
      title(['Grid (W-points) Section at  \xi = ',num2str(index)]);
    end

    if (plt == 2)
      h2=subplot(2,1,2);
      p2=get(h2,'pos');
      p2(4)=p2(4)+0.2;
      set (h2,'pos',p2);

      hold off;
      fill(eta2,hs,[0.6 0.7 0.6]);
      hold on;
      han2=plot(eta',z);
      set(han2,'color', [0.5 0.5 0.5]);
      set(gca,'xlim',[-Inf Inf],'ylim',[zmin 0]);
    end

    xlabel({['Vcoord = ', num2str(G.Vtransform),',',  ...
                          num2str(G.Vstretching),     ...
             '   \theta_s = ' num2str(G.theta_s),     ...
             '   \theta_b = ' num2str(G.theta_b),     ...
             '   hc  = ' num2str(G.hc),               ...
             '   N = ' num2str(G.N)],['y-axis']});
    ylabel('depth  (m)');

  else

    set(gcf,'Units','Normalized',           ...
        'Position',[0.2 0.1 0.6 0.8],       ...
        'PaperOrientation', 'landscape',    ...
        'PaperUnits','Normalized',          ...
        'PaperPosition',[0.2 0.1 0.6 0.8]);

    if (G.spherical)
      xi=G.lon_rho(:,index)';
    else
      xi=G.x_rho(:,index)';
    end
    xi2=[xi(1) xi xi(Lp)];
    hs=-G.h(:,index)';
    zmin=min(hs);
    hs=[zmin hs zmin];

    if (plt == 2)
      h1=subplot(2,1,1);
      p1=get(h1,'pos');
      p1(2)=p1(2)+0.1;
      p1(4)=p1(4)-0.1;
      set (h1,'pos',p1);
    end

    hold off;
    fill(xi2,hs,[0.6 0.7 0.6]);
    hold on;
    han1=plot(xi,z);
    set(han1,'color', [0.5 0.5 0.5]);
    if (plt == 2)
      set(gca,'xlim',[-Inf Inf],'ylim',[-abs(Zzoom) 0]);
      ylabel('depth  (m)');
    else
      set(gca,'xlim',[-Inf Inf],'ylim',[zmin 0]);
    end

    if (kgrid == 0)
      title(['Grid (\rho-points) Section at  \eta = ',num2str(index)]);
    else
      title(['Grid (W-points) Section at  \eta = ',num2str(index)]);
    end

    if (plt == 2)
      h2=subplot(2,1,2);
      p2=get(h2,'pos');
      p2(4)=p2(4)+0.2;
      set (h2,'pos',p2);

      hold off;
      fill(xi2,hs,[0.6 0.7 0.6]);
      hold on;
      han2=plot(xi,z);
      set(han2,'color', [0.5 0.5 0.5]);
      set(gca,'xlim',[-Inf Inf],'ylim',[zmin 0]);
    end

    xlabel({['Vcoord = ', num2str(G.Vtransform),',',  ...
                          num2str(G.Vstretching),     ...
             '   \theta_s = ' num2str(G.theta_s),     ...
             '   \theta_b = ' num2str(G.theta_b),     ...
             '   hc  = ' num2str(G.hc),               ...
             '   N = ' num2str(G.N)],['x-axis']});
    ylabel('depth  (m)');

  end

end

return


