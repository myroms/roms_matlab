function seagrid2roms (theSeagridFile, theRomsFile, theGridTitle)

% SEAGRID2ROMS: Creates ROMS Grid NetCDF file from SeaGrid file.
%
% seagrid2roms('theSeagridFile', 'theRomsFile', 'theGridTitle')
%
% This function creates ROMS Grid NetCDF file, based on the given
% SeaGrid file.  The get/put dialog boxes are invoked where file
% names are absent or given by a wildcard.
%
% On Input:
%
%   theSeagridFile     Seagrid Matlab file name (character string)
%   theRomsFile        ROMS Grid NetCDF file name (character string)
%   theGridTile        Grid title (character string)
%
% Calls:               c_grid (LP, MP, theRomsFile, NewFile)
%                      [umask, vmask, pmask] = uvp_masks(rmask)
%                      nc_write (...)
  
% svn $Id$
%===========================================================================%
%  Copyright (c) 2002-2013 The ROMS/TOMS Group                              %
%    Licensed under a MIT/X style license                                   %
%    See License_ROMS.txt                           Hernan G. Arango        %
%===========================================================================%

EARTH_RADIUS_METERS = 6371.315*1000;   % Equatorial radius.

if (nargin < 1),
  theSeagridFile = '*.mat';
end

if (nargin < 2), 
  theRomsFile = 'roms_grd.nc';
end

if (nargin < 3),
  theGridTitle = char(zeros(1, 128)+abs(' '));
end

if isempty(theSeagridFile) || any(theSeagridFile == '*')
  [f, p] = uigetfile(theSeagridFile, 'Select SeaGrid File:');
  if ~any(f),
    return,
  end
  if (p(end) ~= filesep),
    p(end+1) = filesep;
  end
  theSeagridFile = [p f]
end

if nargin < 2 | isempty(theSeagridFile) | any(theRomsFile == '*')
  [f, p] = uiputfile(theRomsFile, 'Save to Roms File:');
  if ~any(f),
    return,
  end
  if (p(end) ~= filesep),
    p(end+1) = filesep;
  end
  theRomsFile = [p f]
end

disp(' ')
disp([' ## SeaGrid Source File  : ' theSeagridFile])
disp([' ## ROMS Destination File: ' theRomsFile])

% Load the SeaGrid file and get parameters.

try
  theSeagridData = load(theSeagridFile, 's');
catch
  disp([' ## Unable to load: "' theSeagridFile '"'])
  return
end

% With grid_x of size [m, n], the grid itself has [m-1, n-1] cells.
% The latter size corresponds to the size of the mask and bathymetry.
% These cell-centers are called the "rho" points.

s = theSeagridData.s;

grid_x = s.grids{1} * EARTH_RADIUS_METERS;
grid_y = s.grids{2} * EARTH_RADIUS_METERS;
[m, n] = size(grid_x);

geogrid_lon = s.geographic_grids{1};
geogrid_lat = s.geographic_grids{2};
geometry = s.geometry;

% Initialize land/sea masking

mask = s.mask;   % land = 1; water = 0.

if ~isequal(size(mask), size(grid_x)-1)
  if ~isempty(mask)
    disp(' ## Wrong size mask.')
  end
  mask = zeros(m-1, n-1);
end

mask = ~~mask;
land = mask;
water = ~land;

% Load bathymetry

bathymetry = s.gridded_bathymetry;
projection = s.projection;
ang = s.orientation * pi / 180;     % ROMS needs radians
min_depth = s.clipping_depths(1);
max_depth = s.clipping_depths(2);

% Clip Bathymetry

bathymetry(find(isnan(bathymetry))) = min_depth;
bathymetry(bathymetry<min_depth) = min_depth;
bathemetry(bathymetry>max_depth) = max_depth;

spaced_x = s.spaced_grids{1};
spaced_y = s.spaced_grids{2};

% Double the grid-size before proceeding.
%  The grid-cell-centers are termed the "rho" points.

theInterpFcn = 'interp2';
theInterpMethod = 'spline';

grid_x = feval(theInterpFcn, grid_x, 1, theInterpMethod);
grid_y = feval(theInterpFcn, grid_y, 1, theInterpMethod);
geogrid_lon = feval(theInterpFcn, geogrid_lon, 1, theInterpMethod);
geogrid_lat = feval(theInterpFcn, geogrid_lat, 1, theInterpMethod);
spaced_x = feval(theInterpFcn, spaced_x, 1, theInterpMethod);
spaced_y = feval(theInterpFcn, spaced_y, 1, theInterpMethod);

% The present size of the grid nodes.

[n, m] = size(grid_x);

% Flip arrays top for bottom.

FLIPPING = 1;

if FLIPPING
  grid_x = flipud(grid_x);
  grid_y = flipud(grid_y);
  geogrid_lon = flipud(geogrid_lon);
  geogrid_lat = flipud(geogrid_lat);
  geometry{1} = flipud(geometry{1});
  geometry{2} = flipud(geometry{2});
  mask = flipud(mask);
  bathymetry = flipud(bathymetry);
  ang = flipud(ang);
  spaced_x = flipud(spaced_x);
  spaced_y = flipud(spaced_y);
end

xl = max(grid_x(:)) - min(grid_x(:));
el = max(grid_y(:)) - min(grid_y(:));

%---------------------------------------------------------------------------
%  Create ROMS Grid NetCDF file
%---------------------------------------------------------------------------

% The SeaGrid is now a full array, whose height and width are odd-valued.
% We extract staggered sub-grids for the ROMS C-grid type, ignoring the
% outermost rows and columns.  Thus, the so-called "rho" points correspond
% to the even-numbered points in an (i, j) Matlab array.  The "psi" points
% begin at i = 3 and j = 3.  The whole set is indexed as follows:
%
%    rho-points (2:2:end-1, 2:2:end-1), i.e. (2:2:m, 2:2:n)
%    psi-points (3:2:end-2, 3:2:end-2)
%    u-points   (2:2:end-1, 3:2:end-2)
%    v-points   (3:2:end-2, 2:2:end-1)

if ~rem(m, 2),
  m = m-1;        % m must be odd.
end
if ~rem(n, 2),
  n = n-1;        % n must be odd
end

i_rho = 2:2:m-1; j_rho = 2:2:n-1;
i_psi = 3:2:m-2; j_psi = 3:2:n-2;
i_u   = 3:2:m-2; j_u   = 2:2:n-1;
i_v   = 2:2:m-1; j_v   = 3:2:n-2;

% The xi direction (left-right):

LP = (m-1)/2;   % The rho dimension.
L = LP-1;       % The psi dimension.

% The eta direction (up-down):

MP = (n-1)/2;   % The rho dimension.
M = MP-1;       % The psi dimension.

% Create ROMS Grid NetCDF file

NewFile = true;

status = c_grid (LP, MP, theRomsFile, NewFile);

%---------------------------------------------------------------------------
% Extract and set ROMS Grid variables.
%---------------------------------------------------------------------------

switch lower(projection)
  case 'mercator'
    theProjection = 'ME';
  case 'stereographic'
    theProjection = 'ST';
  case 'lambert conformal conic'
    theProjection = 'LC';
end

% Set ROMS Grid variables.

x_rho   = grid_x(j_rho, i_rho);
y_rho   = grid_y(j_rho, i_rho);

x_psi   = grid_x(j_psi, i_psi);
y_psi   = grid_y(j_psi, i_psi);

x_u     = grid_x(j_u, i_u);
y_u     = grid_y(j_u, i_u);

x_v     = grid_x(j_v, i_v);
y_v     = grid_y(j_v, i_v);

lon_rho = geogrid_lon(j_rho, i_rho);
lat_rho = geogrid_lat(j_rho, i_rho);

lon_psi = geogrid_lon(j_psi, i_psi);
lat_psi = geogrid_lat(j_psi, i_psi);

lon_u   = geogrid_lon(j_u, i_u);
lat_u   = geogrid_lat(j_u, i_u);

lon_v   = geogrid_lon(j_v, i_v);
lat_v   = geogrid_lat(j_v, i_v);

% Compute Coriolis parameter.

deg2rad = pi / 180.0;
omega   = 2.0 * pi * 366.25 / (24.0 * 3600.0 * 365.25);

f = 2.0 * omega * sin(deg2rad * lat_rho);

% Set bathymetry.

if (~isempty(bathymetry)),
  h = bathymetry;
end

% Set land/sea masking arrays.

mask = ~~mask;
land = mask;
water = ~land;

rmask = double(water);
rmask = rmask';          % need to transpose here

[umask, vmask, pmask] = uvp_masks(rmask);

% Set curvininear rotation angle.

temp = 0.5*(ang(1:end-1, :) + ang(2:end, :));
ang = zeros(n, m);
ang(2:2:end, 2:2:end) = temp;

angle = ang(j_rho, i_rho);

% Set grid spacing metrics.

if (0)
  sx = abs(spaced_x(:, 2:end) - spaced_x(:, 1:end-1));
  sy = abs(spaced_x(2:end, :) - spaced_x(1:end-1, :));
elseif (0)
  sx = abs(spaced_x(1:2:end, 3:2:end) - spaced_x(1:2:end, 1:2:end-2));
  sy = abs(spaced_y(3:2:end, 1:2:end) - spaced_y(1:2:end-2, 1:2:end));
  sx = 0.5 * (sx(1:end-1, :) + sx(2:end, :));
  sy = 0.5 * (sy(:, 1:end-1) + sy(:, 2:end));
end

% Use geometry from seagrid file.
% Note: need half the number of points.

gx = geometry{1};   % Spherical distances in meters.
gy = geometry{2};

% raw_grid_size  = [m, n]
% geometry_sizes = [size(gx) size(gy)]

sx = 0.5*(gx(1:end-1, :) + gx(2:end, :));
sy = 0.5*(gy(:, 1:end-1) + gy(:, 2:end));

% raw_s_sizes = [size(sx) size(sy)]

% sx = sx(2:end-1, :);
% sy = sy(:, 2:end-1);

pm = 1.0 ./ sx;
pn = 1.0 ./ sy;

% Set curvilinear coordinates inverse metric derivatives.

dmde = zeros(size(pm));
dndx = zeros(size(pn));

dmde(2:end-1, :) = 0.5*(1./pm(3:end, :) - 1./pm(1:end-2, :));
dndx(:, 2:end-1) = 0.5*(1./pn(:, 3:end) - 1./pn(:, 1:end-2));

%---------------------------------------------------------------------------
% Write out ROMS Grid variables into NetCDF file.  We need to transpose
% the data since C-grid computes everything in row-major order (like
% that of the C-language).  We need the arrays in column-major order
% (like Fortran) when using nc_write.
%---------------------------------------------------------------------------

status = nc_write (theRomsFile, 'spherical', 1);


status = nc_write (theRomsFile, 'xl', xl);
status = nc_write (theRomsFile, 'el', el);

status = nc_write (theRomsFile, 'angle', angle');

status = nc_write (theRomsFile, 'pn', pn');
status = nc_write (theRomsFile, 'pm', pm');

status = nc_write (theRomsFile, 'dmde', dmde');
status = nc_write (theRomsFile, 'dndx', dndx');

status = nc_write (theRomsFile, 'f', f');

if (~isempty(bathymetry)),
  rec = 1;
  status = nc_write (theRomsFile, 'hraw', h', rec);
  status = nc_write (theRomsFile, 'h', h');
end
  
status = nc_write (theRomsFile, 'x_rho', x_rho');
status = nc_write (theRomsFile, 'y_rho', y_rho');
status = nc_write (theRomsFile, 'x_psi', x_psi');
status = nc_write (theRomsFile, 'y_psi', y_psi');
status = nc_write (theRomsFile, 'x_u'  , x_u');
status = nc_write (theRomsFile, 'y_u'  , y_u');
status = nc_write (theRomsFile, 'x_v'  , x_v');
status = nc_write (theRomsFile, 'y_v'  , y_v');

status = nc_write (theRomsFile, 'lon_rho', lon_rho');
status = nc_write (theRomsFile, 'lat_rho', lat_rho');
status = nc_write (theRomsFile, 'lon_psi', lon_psi');
status = nc_write (theRomsFile, 'lat_psi', lat_psi');
status = nc_write (theRomsFile, 'lon_u'  , lon_u');
status = nc_write (theRomsFile, 'lat_u'  , lat_u');
status = nc_write (theRomsFile, 'lon_v'  , lon_v');
status = nc_write (theRomsFile, 'lat_v'  , lat_v');

status = nc_write (theRomsFile, 'mask_rho', rmask);  % already transposed
status = nc_write (theRomsFile, 'mask_psi', pmask);  % already transposed
status = nc_write (theRomsFile, 'mask_u',   umask);  % already transposed
status = nc_write (theRomsFile, 'mask_v',   vmask);  % already transposed

return
