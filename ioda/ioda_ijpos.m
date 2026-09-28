function [X,Y]=ioda_ijpos(G, InpName, OutName, varargin)

%
% IODA_IJPOS:  Computes X- and Y-locations in fractional coordinates
%
% [X,Y]=ioda_ijpos(G, InpName, OutName, Lplot, Correction, obc_edge,
%                  Ioffset, Joffset)
%
% This function calculates IODA 'MetaData' observation fractional grid
% locations (x_grid, y_grid) in terms of ROMS (I, J) coordinates. The
% purpose of this computation is to facilitate the processing of H(x)
% operators in ROMS for curvilinear applications.
%
% If the Ioffset and Joffset vectors are provided, the polygon defined
% by the application grid will be reduced by the number of grid points
% specified in the offsets. This adjustment is made to avoid processing
% near the application boundary.
%
% The (x_grid, y_grid) are only used in native ROMS and ignored in
% ROMS-JEDI.  This function is usefull in Nested 4D-Var applications.
%
% On Input:
%
%    G             Input ROMD Grid NetCDF file name (string)
%              or, an existing grid structure (struct)
%
%    InpName       Input IODA NetCDF-4 filename (string)
%
%    OutName       Output thinned IODA NetCDF-4 filename (string)
%
%    Lplot         Switch to plot observation locations (default=false)
%  
%    Correction    Switch to apply correction due to spherical/curvilinear
%                    grids (default=false)
%
%    obc_edge      Switch to include observations on open boundary edges
%                    (default=false)
%
%    Ioffset       Application I-grid offset when defining polygon
%                  (vector; default=[0,0]):
%                    Ioffset(1):  I-grid offset on the edge where Istr=1
%                    Ioffset(2):  I-grid offset on the edge where Iend=Lm
%
%    Joffset       Application J-grid offset when defining polygon
%                  (vector; default=[0,0]):
%                    Joffset(1):  J-grid offset on the edge where Jstr=1
%                    Joffset(2):  J-grid offset on the edge where Jend=Mm
%
% On Output:
%
%    X             Observation fractional x-grid location
%
%    Y             Observation fractional x-grid location
%

% git $Id$
%=======================================================================%
%  Copyright (c) 2002-2026 The ROMS Group                               %
%    Licensed under a MIT/X style license                               %
%    See License_ROMS.md                            Hernan G. Arango    %
%=======================================================================%

switch numel(varargin)
  case 0
   Lplot = false;
   Correction = false;
    obc_edge = false;
    Ioffset = [0, 0];
    Joffset = [0, 0];
  case 1
    Lplot = varargin{1};
    Correction = false;
    obc_edge = false;
    Ioffset = [0, 0];
    Joffset = [0, 0];
  case 2
    Lplot = varargin{1};
    Correction = varargin{2};
    obc_edge = false;
    Ioffset = [0, 0];
    Joffset = [0, 0];
  case 3
    Lplot = varargin{1};
    Correction = varargin{2};
    obc_edge = varargin{3};
    Ioffset = [0, 0];
    Joffset = [0, 0];
  case 4
    Lplot = varargin{1};
    Correction = varargin{2};
    obc_edge = varargin{3};
    Ioffset = varargin{4};
    Joffset = [0, 0];
  case 5
    Lplot = varargin{1};
    Correction = varargin{2};
    obc_edge = varargin{3};
    Ioffset = varargin{4};
    Joffset = varargin{5};
end

if (ischar(G))
  G = get_roms_grid(G);
end

%------------------------------------------------------------------------
% Read in input IODA NetCDF-4 observation file.
%------------------------------------------------------------------------

Obs = ioda_read(InpName);

%------------------------------------------------------------------------
%  Extract polygon defining application grid box.
%------------------------------------------------------------------------

%  Set grid application polygon.

[Im,Jm]=size(G.lon_rho);

Istr = 1 +Ioffset(1);
Iend = Im-Ioffset(2);
Jstr = 1 +Joffset(1);
Jend = Jm-Joffset(2);

Xbox = [squeeze(G.lon_rho(Istr:Iend,Jstr));                           ...
        squeeze(G.lon_rho(Iend,Jstr+1:Jend))';                        ...
        squeeze(flipud(G.lon_rho(Istr:Iend-1,Jend)));                 ...
        squeeze(fliplr(G.lon_rho(Istr,Jstr:Jend-1)))'];

Ybox = [squeeze(G.lat_rho(Istr:Iend,Jstr));                           ...
        squeeze(G.lat_rho(Iend,Jstr+1:Jend))';                        ...
        squeeze(flipud(G.lat_rho(Istr:Iend-1,Jend)));                 ...
        squeeze(fliplr(G.lat_rho(Istr,Jstr:Jend-1)))'];

%  Find observation inside (IN) or on the edge (ON) the polygon defined
%  by (Xbox,Ybox).

[IN ON] = inpolygon(Obs.longitude, Obs.latitude, Xbox, Ybox);

%  Flag outlier observations as bounded=false.  We are only considering
%  observations inside the polygon.

bounded=false(size(Obs.longitude));

bounded(IN) = true;
if (obc_edge)                  % process observations on boundary edges
  bounded(ON) = true;
end

%------------------------------------------------------------------------
% Plot observations in the application grid.
%------------------------------------------------------------------------

if (Lplot)
  pcolor(G.lon_rho, G.lat_rho, ones(size(G.lon_rho)));
  hold on;
  set(gca, 'fontsize', 14, 'fontweight', 'bold');
  if (isfield(G, 'lon_coast'))
    plot(G.lon_coast, G.lat_coast, 'k-');
  end
  plot(Xbox, Ybox, 'b-',                                              ...
       Obs.longitude(IN), Obs.latitude(IN), 'b.',                     ...
       Obs.longitude(ON), Obs.latitude(ON), 'r.',                     ...
       Obs.longitude(~IN), Obs.latitude(~IN), 'm.');
  title(['Blue points (inside),  ',                                   ...
         'Red points (boundary),  ',                                  ...
         'Magenta Points (outliers)'],                                ...
         'fontsize', 14, 'fontweight', 'bold');
  hold off;
end

%------------------------------------------------------------------------
% Compute model grid fractional (I,J) locations at observation locations
% via interpolation.
%------------------------------------------------------------------------

Igrid =repmat([0:1:Im-1]', [1 Jm]);
Jgrid =repmat([0:1:Jm-1] , [Im 1]);

if (isfield(G, 'mask_rho'))
  ind = find(G.mask_rho < 0.5);
  if (~isempty(ind))
    Igrid(ind) = NaN;
    Jgrid(ind) = NaN;
  end
end

% Initialize unbounded observations to NaN.

X = ones(size(Obs.longitude)) .* NaN;
Y = ones(size(Obs.latitude)) .* NaN;

X(bounded) = griddata(G.lon_rho, G.lat_rho, Igrid,                    ...
                      Obs.longitude(bounded), Obs.latitude(bounded));
Y(bounded) = griddata(G.lon_rho, G.lat_rho, Jgrid,                    ...
                      Obs.longitude(bounded), Obs.latitude(bounded));

%  If land/sea masking, find the observation in land (Xgrid=Ygrid=NaN);

ind = find(isnan(X) | isnan(Y));
if (~isempty(ind))
  bounded(ind) = false;
end

%------------------------------------------------------------------------
%  Spherical/Curvilinear corrections.
%------------------------------------------------------------------------

if (Correction)
  [X, Y] = correction(G.lon_rho, G.lat_rho, G.angle,                  ...
                      Obs.longitude, Obs.latitude, bounded, X, Y);
end

%------------------------------------------------------------------------
%  Create and write output IODA NetCDF-4 file.
%------------------------------------------------------------------------

roms_LmMmN = Obs.roms_grid;
roms_LmMmN(1) = G.Lm;
roms_LmMmN(2) = G.Mm;
Obs.roms_grid = roms_LmMmN;        % in case of nested grid observations

Obs.ncname = OutName;
Obs.x_grid = X;
Obs.y_grid = Y;

create_ioda_obs(Obs, OutName);
ioda_write(Obs, OutName);

return


function [X,Y]=correction(rlon, rlat, angle, obs_lon, obs_lat,        ...
                          bounded, X, Y)

%
% CORRECTION:  Apply curvilinear coordinates correction to observations
%              fractional (X,Y) locations.
%
% The maximum variable size allowed by Matlab can be exceeded very quickly
% if the observation vector is large. This function is coded either using
% block temporary arrays or a simple loop where the corrections are computed
% one by one. This is one of the few instances in Matlab that actually is
% more efficient to do complex scalar operations than vector operations.
% Notice that the number of satellite observation can be large and often
% exceeds 1E5.
%
% On Input:
%
%    rlon          Application grid longitude at RHO-points (matrix)
%    rlat          Application grid latitude  at RHO-points (matrix)
%    angle         curvilinear grid rotation (radians; matrix)
%    obs_lon       Observation longitude locations (vector)
%    obs_lat       Observation latitude  locations (vector)
%    bounded       Switch marking outlier (false) points (vector)
%    X             Observation fractional x-grid location (first guess)
%    Y             Observation fractional y-grid location (first guess)
%
% On Output:
%
%    X             Adjusted observation fractional x-grid location
%    Y             Adjusted observation fractional y-grid location
%

%=======================================================================%
%  Copyright (c) 2002-2026 The ROMS Group                               %
%    Licensed under a MIT/X style license                               %
%    See License_ROMS.md                            Hernan G. Arango    %
%=======================================================================%

debugging = false;                  % debugging switch

Nobs = length(obs_lon);             % number of observations

block_length = 100000;              % size of block arrays

Eradius = 6371315.0;                % Earth radius (meters)
deg2rad = pi/180;                   % degrees to radians factor

[Lr, Mr] = size(rlon);              % Number of rho-points

Lm = Lr - 1;
Mm = Mr - 1;

%------------------------------------------------------------------------
%  Set (I,J) coordinates of the grid cell containing the observation
%  need to add 1 because zero lower bound in ROMS "rlon"
%------------------------------------------------------------------------

I = fix(X);
ind = find((0 <= I) & (I < Lm));
if (~isempty(ind))
  I(ind) = I(ind) + 1;
end

J = fix(Y);
ind = find((0 <= J) & (J < Mm));
if (~isempty(ind))
  J(ind) = J(ind) + 1;
end

if (debugging)
  disp(' ');
  disp(['  Xmin = ', num2str(min(X),'%6.2f'),                         ...
        '  Xmax = ', num2str(max(X),'%6.2f'),                         ...
        '  Imin = ', num2str(min(I),'%3.3i'),                         ...
        '  Imax = ', num2str(max(I),'%3.3i'),                         ...
        '  Lr = ',   num2str(Lr,'%3.3i')]);

  disp(['  Ymin = ', num2str(min(Y),'%6.2f'),                         ...
        '  Ymax = ', num2str(max(Y),'%6.2f'),                         ...
        '  Jmin = ', num2str(min(J),'%3.3i'),                         ...
        '  Jmax = ', num2str(max(J),'%3.3i'),                         ...
        '  Mr = ',   num2str(Mr,'%3.3i')]);
  disp(' ');
end

%------------------------------------------------------------------------
%  It is possible that we are processing a large number of observations.
%  Therefore, the observation vector is processed by blocks to reduce
%  the memory requirements.
%------------------------------------------------------------------------

N  = ceil(Nobs / block_length);
n1 = 0;
n2 = 0;

while (n2 < Nobs)

  n1 = n2 + 1;
  n2 = n1 + block_length;
  if (n2 > Nobs)
    n2 = Nobs;
  end

  iobs = n1:n2;
  ind  = find(~bounded(n1:n2));
  if (~isempty(ind))              % remove unbouded observations, if any
    iobs(ind) = [];
  end

  i_j   = sub2ind(size(rlon), I(iobs)  , J(iobs)  );
  ip1_j = sub2ind(size(rlon), I(iobs)+1, J(iobs)  );
  i_jp1 = sub2ind(size(rlon), I(iobs)  , J(iobs)+1);

  if (debugging)
    disp(['  Processing observation vector, n1:n2 = '                 ...
          num2str(n1,'%7.7i'), ' - ', num2str(n2,'%7.7i')             ...
          '  size = ', num2str(length(iobs))]);
  end

%  Convert all positions to meters first.

  yfac = Eradius * deg2rad;
  xfac = yfac .* cos(obs_lat(iobs) .* deg2rad);

  xpp  = (obs_lon(iobs) - rlon(i_j)) .* xfac;
  ypp  = (obs_lat(iobs) - rlat(i_j)) .* yfac;

%  Use Law of Cosines to get cell parallelogram "shear" angle.

  diag2 = (rlon(ip1_j) - rlon(i_jp1)) .^ 2 +                          ...
          (rlat(ip1_j) - rlat(i_jp1)) .^ 2;

  aa2   = (rlon(i_j)   - rlon(ip1_j)) .^ 2 +                          ...
          (rlat(i_j)   - rlat(ip1_j)) .^ 2;

  bb2   = (rlon(i_j)   - rlon(i_jp1)) .^ 2 +                          ...
          (rlat(i_j)   - rlat(i_jp1)) .^ 2;

  phi   = asin((diag2 - aa2- bb2) ./ (2 .* sqrt(aa2 .* bb2)));

%  Transform fractional locations into curvilinear coordinates. Assume
%  the cell is rectanglar, for now.

  dx = xpp .* cos(angle(i_j)) + ypp .* sin(angle(i_j));

  dy = ypp .* cos(angle(i_j)) - xpp .* sin(angle(i_j));

%  Correct for parallelogram.

  dx = dx + dy .* tan(phi);
  dy = dy ./ cos(phi);

%  Scale with cell side lengths to translate into cell indexes.

  dx = dx ./ sqrt(aa2) ./ xfac;

  ind = find(dx < 0);
  if (~isempty(ind))
    dx(ind) = 0;
  end

  ind = find(dx > 1);
  if (~isempty(ind))
    dx(ind) = 1;
  end

  dy = dy ./ sqrt(bb2) ./ yfac;

  ind = find(dy < 0);
  if (~isempty(ind))
    dy(ind) = 0;
  end

  ind = find(dy > 1);
  if (~isempty(ind))
    dy(ind) = 1;
  end

  X(iobs) = fix(X(iobs)) + dx;
  Y(iobs) = fix(Y(iobs)) + dy;

end

return