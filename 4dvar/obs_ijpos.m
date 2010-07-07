function [Xgrid,Ygrid]=obs_ijpos(GRDname,obs_lon,obs_lat,Ccorrection);

%
% OBS_IJPOS:  Computes observation locations in ROMS fractional coordinates
%
% [Xgrid,Ygrid]=obs_ijpos(GRDname,obs_lon,obs_lat,Ccorrection)
%
% This function computes the observation locations (Xgrid,Ygrid) in terms
% of ROMS fractional (I,J) coordinates. This is done to facilitate the
% processing of the observation operators inside ROMS. All the observations
% are assumed to be located at RHO-points.
%
% On Input:
%
%    GRDname       NetCDF grid file name (string)
%    obs_lon       observation longitude (positive, degrees_east)
%    obs_lat       observation latitude  (positive, degrees_north)
%    Ccorrection   switch to apply small correction in curvilinear
%                    grids (0: no, 1: yes)
%
% On Ouput:
%
%    Xgrid         observation fractional x-grid location
%    Ygrid         observation fractional y-grid location
%
% Notice:  Outlier observations outside of ROMS grid has an NaN value in
%          Xgrid and Ygrid.
%
  
% svn $Id$
%===========================================================================%
%  Copyright (c) 2002-2010 The ROMS/TOMS Group                              %
%    Licensed under a MIT/X style license                                   %
%    See License_ROMS.txt                           Hernan G. Arango        %
%===========================================================================%

IPLOT=0;          % switch for plotting during debugging

% If appropriate, transpose input vectors so they are of size: [length 1].

lsize=size(obs_lon);
if (lsize(1) < lsize(2)),
  obs_lon=obs_lon';
end,

lsize=size(obs_lat);
if (lsize(1) < lsize(2)),
  obs_lat=obs_lat';
end,

if (nargin < 4),
  Ccorrection=0;
end,

%----------------------------------------------------------------------------
% Read in model grid variables.
%----------------------------------------------------------------------------

[vname,nvars]=nc_vname(GRDname);

got.spherical=0;
got.angle=0;
got.mask_rho=0;
got.lon_rho=0;
got.lat_rho=0;

for n=1:nvars
  name=deblank(vname(n,:));
  switch name
    case 'spherical'
      got.spherical=1;
    case 'lon_rho'
      got.lon_rho=1;
    case 'lat_rho'
      got.lat_rho=1;
    case 'angle'
      got.angle=1;
    case 'mask_rho'
      got.mask_rho=1;
  end,
end,

if (got.spherical),
  spherical=nc_read(GRDname,'spherical');
  if (ischar(spherical)),
    if (spherical == 'T' | spherical == 't'),
      spherical=1;
    else,
      spherical=0;
    end,
  end,
else,
  spherical=1;
end,

if (got.lon_rho),
  rlon=nc_read(GRDname,'lon_rho');
else,
  error(['OBS_IJPOS - cannot find variable: lon_rho']);
end,

if (got.lat_rho),
  rlat=nc_read(GRDname,'lat_rho');
else,
  error(['OBS_IJPOS - cannot find variable: lat_rho']);
end,

if (got.angle),
  angle=nc_read(GRDname,'angle');
else
  angle=zeros(size(rlon));
end,

if (got.mask_rho),
  rmask=nc_read(GRDname,'mask_rho');
else
  rmask=ones(size(rlon));
end,

%----------------------------------------------------------------------------
%  Extract polygon defining application grid box.
%----------------------------------------------------------------------------

[Im,Jm]=size(rlon);

Xbox=[squeeze(rlon(:,1)); ...
      squeeze(rlon(Im,2:Jm))'; ...
      squeeze(flipud(rlon(1:Im-1,Jm))); ...
      squeeze(fliplr(rlon(1,1:Jm-1)))'];

Ybox=[squeeze(rlat(:,1)); ...
      squeeze(rlat(Im,2:Jm))'; ...
      squeeze(flipud(rlat(1:Im-1,Jm))); ...
      squeeze(fliplr(rlat(1,1:Jm-1)))'];

%  Set complex vector of points outlining polygon.

Zbox=complex(Xbox,Ybox);

%  Check if requested observation locations are inside/outside of the
%  outlining polygon. Outliers are marked as k=NaN.  The 'inside' routine
%  returns also an NaN if the observation is right on the boundary
%  (Xbox,Ybox).  Perhaps, we need to find a way to fix this routine.

p=complex(obs_lon,obs_lat); 
k=inside(p,Zbox);                   % gives indexes of inside points only

outlier=ones(size(obs_lon)).*NaN;   % flag outside indexes with NaN
if (~isempty(k)),
  for n=1:length(k),
    m=k(n);
    outlier(m)=m;
  end,
end,

% Plot observations in the application grid.

if (IPLOT),
  pcolor(rlon,rlat,ones(size(rlon)));
  hold on;
  plot(obs_lon,obs_lat,'b.');
  hold off;
end,

%----------------------------------------------------------------------------
% Compute model grid fractional (I,J) locations at observation locations
% via interpolation.
%----------------------------------------------------------------------------

Igrid=repmat([0:1:Im-1]',[1 Jm]);
Jgrid=repmat([0:1:Jm-1] ,[Im 1]);

if (got.mask_rho),
  ind=find(rmask < 1);
  if (~isempty(ind)),
    Igrid(ind)=NaN;
    Jgrid(ind)=NaN;
  end,
end,

Xgrid=griddata(rlon,rlat,Igrid,obs_lon,obs_lat);
Ygrid=griddata(rlon,rlat,Jgrid,obs_lon,obs_lat);

%  If land/sea masking, find the observation in land (Xgrid=Ygrid=NaN);

if (got.mask_rho),
  ind=find(isnan(Xgrid) | isnan(Ygrid));
  if (~isempty(ind)),
    outlier(ind)=NaN;
  end,
end,
outlier=find(isnan(outlier));

%----------------------------------------------------------------------------
% Curvilinear corrections.
%----------------------------------------------------------------------------

if (Ccorrection),

  Eradius=6371315.0;                % meters
  deg2rad=pi/180;                   % degrees to radians factor

  I=fix(Xgrid)+1;                   % need to add 1 because zero lower bound
  J=fix(Ygrid)+1;                   % need to add 1 because zero lower bound

% If outliers, set (I,J) indexes to (1,1) for now to allow the
% vector computations below. Their values will be marked with
% NaN at the end. The curvilinear correction is irrelevant in
% those cases.

if (~isempty(outlier));
  I(outlier)=1;
  J(outlier)=1;
end,
  
% Knowing the correct cell, calculate the exact indexes, accounting
% for a possibly rotated grid.  Convert all positions to meters first.

  yfac=Eradius*deg2rad;
  xfac=yfac.*cos(obs_lat.*deg2rad);

  xpp=(obs_lon-diag(rlon(I,J))).*xfac;
  ypp=(obs_lat-diag(rlat(I,J))).*yfac;

% Use Law of Cosines to get cell parallelogram "shear" angle.

  diag2=(diag(rlon(I+1,J))-diag(rlon(I,J+1))).^2+ ...
        (diag(rlat(I+1,J))-diag(rlat(I,J+1))).^2;

  aa2=(diag(rlon(I,J))-diag(rlon(I+1,J))).^2+ ...
      (diag(rlat(I,J))-diag(rlat(I+1,J))).^2;

  bb2=(diag(rlon(I,J))-diag(rlon(I,J+1))).^2+ ...
      (diag(rlat(I,J))-diag(rlat(I,J+1))).^2;

  phi=asin((diag2-aa2-bb2)./(2.*sqrt(aa2.*bb2)));

% Transform fractional locations into curvilinear coordinates. Assume the
% cell is rectanglar, for now.

  dx=xpp.*cos(diag(angle(I,J))+ypp.*sin(diag(angle(I,J))));
  dy=ypp.*cos(diag(angle(I,J))-xpp.*sin(diag(angle(I,J))));

% Correct for parallelogram.

  dx=dx+dy.*tan(phi);
  dy=dy./cos(phi);

% Scale with cell side lengths to translate into cell indexes.

  dx=dx./sqrt(aa2)./xfac;

  ind=find(dx < 0); if (~isempty(ind)), dx(ind)=0; end,
  ind=find(dx > 1); if (~isempty(ind)), dx(ind)=1; end,

  dy=dy./sqrt(bb2)./yfac;

  ind=find(dy < 0); if (~isempty(ind)), dy(ind)=0; end,
  ind=find(dy > 1); if (~isempty(ind)), dy(ind)=1; end,

  Xgrid=fix(Xgrid)+dx;
  Ygrid=fix(Ygrid)+dy;

end,

%----------------------------------------------------------------------------
% Mark outlier obsvations.
%----------------------------------------------------------------------------

if (~isempty(outlier));
  Xgrid(outlier)=NaN;
  Ygrid(outlier)=NaN;
end,

return
