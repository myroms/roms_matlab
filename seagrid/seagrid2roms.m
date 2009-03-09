function seagrid2roms(theSeagridFile, theRomsFile, theGridTitle)

% seagrid2roms -- Output ROMS format from SeaGrid file.
%  seagrid2roms('theSeagridFile', 'theRomsFile', 'theGridTitle')
%   creates a ROMS file, based on the given SeaGrid file.  The
%   get/put dialog boxes are invoked where filenames are absent
%   or given by a wildcard.
%
% NOTE: this routine does not employ _FillValue attributes
%  in the output NetCDF variables.
 
% svn $Id$
%=======================================================================
% Copyright (C) 1999 Dr. Charles R. Denham, ZYDECO.
%  All Rights Reserved.
%   Disclosure without explicit written consent from the
%    copyright owner does not constitute publication.
%=======================================================================
 
% Version of 18-Jun-1999 17:11:59.
% Updated    23-Oct-2000 11:30:58.

RADIAN_CONVERSION_FACTOR = 180/pi;
EARTH_RADIUS_METERS = 6378*1000;   % Equatorial radius.

if nargin < 1, theSeagridFile = '*.mat'; end
if nargin < 2, theRomsFile = 'roms_grd.nc'; end
if nargin < 3, theGridTitle = char(zeros(1, 128)+abs(' ')); end

if isempty(theSeagridFile) | any(theSeagridFile == '*')
		 [f, p] = uigetfile(theSeagridFile, 'Select SeaGrid File:');
		 if ~any(f), return, end
		 if p(end) ~= filesep, p(end+1) = filesep; end
		 theSeagridFile = [p f]
end

if nargin < 2 | isempty(theSeagridFile) | any(theRomsFile == '*')
		 [f, p] = uiputfile(theRomsFile, 'Save to Roms File:');
		 if ~any(f), return, end
		 if p(end) ~= filesep, p(end+1) = filesep; end
		 theRomsFile = [p f]
end

disp([' ## SeaGrid Source File  : ' theSeagridFile])
disp([' ## ROMS Destination File: ' theRomsFile])

% Load the SeaGrid file and get parameters.

try
		 theSeagridData = load(theSeagridFile, 's');
catch
		 disp([' ## Unable to load: "' theSeagridFile '"'])
		 return
end

% With grid_x of size [m, n], the grid itself has
%  [m-1, n-1] cells.  The latter size corresponds
%  to the size of the mask and bathymetry.  These
%  cell-centers are called the "rho" points.

s = theSeagridData.s;

grid_x = s.grids{1} * EARTH_RADIUS_METERS;
grid_y = s.grids{2} * EARTH_RADIUS_METERS;
[m, n] = size(grid_x);

geogrid_lon = s.geographic_grids{1};
geogrid_lat = s.geographic_grids{2};
geometry = s.geometry;

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

bathymetry = s.gridded_bathymetry;
projection = s.projection;
ang = s.orientation * pi / 180;   % ROMS needs radians.
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
 		 % f = flipud(f);
		 spaced_x = flipud(spaced_x);
		 spaced_y = flipud(spaced_y);
end

xl = max(grid_x(:)) - min(grid_x(:));
el = max(grid_y(:)) - min(grid_y(:));

% Create the Roms NetCDF file.

%% ncdump('priapus:MATLAB:toolbox:Roms:tasman_grd.nc')   %% Generated 17-Jun-1999 17:37:30
 
nc = netcdf(theRomsFile, 'clobber');
if isempty(nc), return, end
 
%% Global attributes:

disp(' ## Defining Global Attributes...')
 
nc.type = ncchar('Gridpak file');
nc.gridid = theGridTitle;
nc.history = ncchar(['Created by "' mfilename '" on ' datestr(now)]);

nc.CPP_options = ncchar('DCOMPLEX, DBLEPREC, NCARG_32, PLOTS,');
name(nc.CPP_options, 'CPP-options')

% The SeaGrid is now a full array, whose height
%  and width are odd-valued.  We extract staggered
%  sub-grids for the Roms scheme, ignoring the
%  outermost rows and columns.  Thus, the so-called
%  "rho" points correspond to the even-numbered points
%  in an (i, j) Matlab array.  The "psi" points begin
%  at i = 3 and j = 3.  The whole set is indexed as
%  follows:

% rho (2:2:end-1, 2:2:end-1), i.e. (2:2:m, 2:2:n), etc.
% psi (3:2:end-2, 3:2:end-2)
% u   (2:2:end-1, 3:2:end-2)
% v   (3:2:end-2, 2:2:end-1)

if ~rem(m, 2), m = m-1; end   % m, n must be odd.
if ~rem(n, 2), n = n-1; end

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

disp(' ## Defining Dimensions...')
 
nc('xi_psi') = L;
nc('xi_rho') = LP;
nc('xi_u') = L;
nc('xi_v') = LP;

nc('eta_psi') = M;
nc('eta_rho') = MP;
nc('eta_u') = MP;
nc('eta_v') = M;

nc('two') = 2;
nc('bath') = 0; %% (record dimension)
 
%% Variables and attributes:

disp(' ## Defining Variables and Attributes...')
 
nc{'xl'} = ncdouble; %% 1 element.
nc{'xl'}.long_name = ncchar('domain length in the XI-direction');
nc{'xl'}.units = ncchar('meter');
 
nc{'el'} = ncdouble; %% 1 element.
nc{'el'}.long_name = ncchar('domain length in the ETA-direction');
nc{'el'}.units = ncchar('meter');
 
nc{'JPRJ'} = ncchar('two'); %% 2 elements.
nc{'JPRJ'}.long_name = ncchar('Map projection type');

nc{'JPRJ'}.option_ME_ = ncchar('Mercator');
nc{'JPRJ'}.option_ST_ = ncchar('Stereographic');
nc{'JPRJ'}.option_LC_ = ncchar('Lambert conformal conic');
name(nc{'JPRJ'}.option_ME_, 'option(ME)')
name(nc{'JPRJ'}.option_ST_, 'option(ST)')
name(nc{'JPRJ'}.option_LC_, 'option(LC)')
 
nc{'PLAT'} = ncfloat('two'); %% 2 elements.
nc{'PLAT'}.long_name = ncchar('Reference latitude(s) for map projection');
nc{'PLAT'}.units = ncchar('degree_north');
 
nc{'PLONG'} = ncfloat; %% 1 element.
nc{'PLONG'}.long_name = ncchar('Reference longitude for map projection');
nc{'PLONG'}.units = ncchar('degree_east');
 
nc{'ROTA'} = ncfloat; %% 1 element.
nc{'ROTA'}.long_name = ncchar('Rotation angle for map projection');
nc{'ROTA'}.units = ncchar('degree');
 
nc{'JLTS'} = ncchar('two'); %% 2 elements.
nc{'JLTS'}.long_name = ncchar('How limits of map are chosen');

nc{'JLTS'}.option_CO_ = ncchar('P1, .. P4 define two opposite corners ');
nc{'JLTS'}.option_MA_ = ncchar('Maximum (whole world)');
nc{'JLTS'}.option_AN_ = ncchar('Angles - P1..P4 define angles to edge of domain');
nc{'JLTS'}.option_LI_ = ncchar('Limits - P1..P4 define limits in u,v space');
name(nc{'JLTS'}.option_CO_, 'option(CO)')
name(nc{'JLTS'}.option_MA_, 'option(MA)')
name(nc{'JLTS'}.option_AN_, 'option(AN)')
name(nc{'JLTS'}.option_LI_, 'option(LI)')
 
nc{'P1'} = ncfloat; %% 1 element.
nc{'P1'}.long_name = ncchar('Map limit parameter number 1');
 
nc{'P2'} = ncfloat; %% 1 element.
nc{'P2'}.long_name = ncchar('Map limit parameter number 2');
 
nc{'P3'} = ncfloat; %% 1 element.
nc{'P3'}.long_name = ncchar('Map limit parameter number 3');
 
nc{'P4'} = ncfloat; %% 1 element.
nc{'P4'}.long_name = ncchar('Map limit parameter number 4');
 
nc{'XOFF'} = ncfloat; %% 1 element.
nc{'XOFF'}.long_name = ncchar('Offset in x direction');
nc{'XOFF'}.units = ncchar('meter');
 
nc{'YOFF'} = ncfloat; %% 1 element.
nc{'YOFF'}.long_name = ncchar('Offset in y direction');
nc{'YOFF'}.units = ncchar('meter');
 
nc{'depthmin'} = ncshort; %% 1 element.
nc{'depthmin'}.long_name = ncchar('Shallow bathymetry clipping depth');
nc{'depthmin'}.units = ncchar('meter');
 
nc{'depthmax'} = ncshort; %% 1 element.
nc{'depthmax'}.long_name = ncchar('Deep bathymetry clipping depth');
nc{'depthmax'}.units = ncchar('meter');
 
nc{'spherical'} = ncchar; %% 1 element.
nc{'spherical'}.long_name = ncchar('Grid type logical switch');
nc{'spherical'}.option_T_ = ncchar('spherical');
nc{'spherical'}.option_F_ = ncchar('Cartesian');
name(nc{'spherical'}.option_T_, 'option(T)')
name(nc{'spherical'}.option_F_, 'option(F)')
 
nc{'hraw'} = ncdouble('bath', 'eta_rho', 'xi_rho'); %% 0 elements.
nc{'hraw'}.long_name = ncchar('Working bathymetry at RHO-points');
nc{'hraw'}.units = ncchar('meter');
nc{'hraw'}.field = ncchar('bath, scalar');
 
nc{'h'} = ncdouble('eta_rho', 'xi_rho'); %% 16900 elements.
nc{'h'}.long_name = ncchar('Final bathymetry at RHO-points');
nc{'h'}.units = ncchar('meter');
nc{'h'}.field = ncchar('bath, scalar');
 
nc{'f'} = ncdouble('eta_rho', 'xi_rho'); %% 16900 elements.
nc{'f'}.long_name = ncchar('Coriolis parameter at RHO-points');
nc{'f'}.units = ncchar('second-1');
nc{'f'}.field = ncchar('Coriolis, scalar');
 
nc{'pm'} = ncdouble('eta_rho', 'xi_rho'); %% 16900 elements.
nc{'pm'}.long_name = ncchar('curvilinear coordinate metric in XI');
nc{'pm'}.units = ncchar('meter-1');
nc{'pm'}.field = ncchar('pm, scalar');
 
nc{'pn'} = ncdouble('eta_rho', 'xi_rho'); %% 16900 elements.
nc{'pn'}.long_name = ncchar('curvilinear coordinate metric in ETA');
nc{'pn'}.units = ncchar('meter-1');
nc{'pn'}.field = ncchar('pn, scalar');
 
nc{'dndx'} = ncdouble('eta_rho', 'xi_rho'); %% 16900 elements.
nc{'dndx'}.long_name = ncchar('xi derivative of inverse metric factor pn');
nc{'dndx'}.units = ncchar('meter');
nc{'dndx'}.field = ncchar('dndx, scalar');
 
nc{'dmde'} = ncdouble('eta_rho', 'xi_rho'); %% 16900 elements.
nc{'dmde'}.long_name = ncchar('eta derivative of inverse metric factor pm');
nc{'dmde'}.units = ncchar('meter');
nc{'dmde'}.field = ncchar('dmde, scalar');
 
nc{'x_rho'} = ncdouble('eta_rho', 'xi_rho'); %% 16900 elements.
nc{'x_rho'}.long_name = ncchar('x location of RHO-points');
nc{'x_rho'}.units = ncchar('meter');
 
nc{'y_rho'} = ncdouble('eta_rho', 'xi_rho'); %% 16900 elements.
nc{'y_rho'}.long_name = ncchar('y location of RHO-points');
nc{'y_rho'}.units = ncchar('meter');
 
nc{'x_psi'} = ncdouble('eta_psi', 'xi_psi'); %% 16641 elements.
nc{'x_psi'}.long_name = ncchar('x location of PSI-points');
nc{'x_psi'}.units = ncchar('meter');
 
nc{'y_psi'} = ncdouble('eta_psi', 'xi_psi'); %% 16641 elements.
nc{'y_psi'}.long_name = ncchar('y location of PSI-points');
nc{'y_psi'}.units = ncchar('meter');
 
nc{'x_u'} = ncdouble('eta_u', 'xi_u'); %% 16770 elements.
nc{'x_u'}.long_name = ncchar('x location of U-points');
nc{'x_u'}.units = ncchar('meter');
 
nc{'y_u'} = ncdouble('eta_u', 'xi_u'); %% 16770 elements.
nc{'y_u'}.long_name = ncchar('y location of U-points');
nc{'y_u'}.units = ncchar('meter');
 
nc{'x_v'} = ncdouble('eta_v', 'xi_v'); %% 16770 elements.
nc{'x_v'}.long_name = ncchar('x location of V-points');
nc{'x_v'}.units = ncchar('meter');
 
nc{'y_v'} = ncdouble('eta_v', 'xi_v'); %% 16770 elements.
nc{'y_v'}.long_name = ncchar('y location of V-points');
nc{'y_v'}.units = ncchar('meter');
 
nc{'lat_rho'} = ncdouble('eta_rho', 'xi_rho'); %% 16900 elements.
nc{'lat_rho'}.long_name = ncchar('latitude of RHO-points');
nc{'lat_rho'}.units = ncchar('degree_north');
 
nc{'lon_rho'} = ncdouble('eta_rho', 'xi_rho'); %% 16900 elements.
nc{'lon_rho'}.long_name = ncchar('longitude of RHO-points');
nc{'lon_rho'}.units = ncchar('degree_east');
 
nc{'lat_psi'} = ncdouble('eta_psi', 'xi_psi'); %% 16641 elements.
nc{'lat_psi'}.long_name = ncchar('latitude of PSI-points');
nc{'lat_psi'}.units = ncchar('degree_north');
 
nc{'lon_psi'} = ncdouble('eta_psi', 'xi_psi'); %% 16641 elements.
nc{'lon_psi'}.long_name = ncchar('longitude of PSI-points');
nc{'lon_psi'}.units = ncchar('degree_east');
 
nc{'lat_u'} = ncdouble('eta_u', 'xi_u'); %% 16770 elements.
nc{'lat_u'}.long_name = ncchar('latitude of U-points');
nc{'lat_u'}.units = ncchar('degree_north');
 
nc{'lon_u'} = ncdouble('eta_u', 'xi_u'); %% 16770 elements.
nc{'lon_u'}.long_name = ncchar('longitude of U-points');
nc{'lon_u'}.units = ncchar('degree_east');
 
nc{'lat_v'} = ncdouble('eta_v', 'xi_v'); %% 16770 elements.
nc{'lat_v'}.long_name = ncchar('latitude of V-points');
nc{'lat_v'}.units = ncchar('degree_north');
 
nc{'lon_v'} = ncdouble('eta_v', 'xi_v'); %% 16770 elements.
nc{'lon_v'}.long_name = ncchar('longitude of V-points');
nc{'lon_v'}.units = ncchar('degree_east');
 
nc{'mask_rho'} = ncdouble('eta_rho', 'xi_rho'); %% 16900 elements.
nc{'mask_rho'}.long_name = ncchar('mask on RHO-points');
nc{'mask_rho'}.option_0_ = ncchar('land');
nc{'mask_rho'}.option_1_ = ncchar('water');
name(nc{'mask_rho'}.option_0_, 'option(0)')
name(nc{'mask_rho'}.option_1_, 'option(1)')
 
nc{'mask_u'} = ncdouble('eta_u', 'xi_u'); %% 16770 elements.
nc{'mask_u'}.long_name = ncchar('mask on U-points');
nc{'mask_u'}.option_0_ = ncchar('land');
nc{'mask_u'}.option_1_ = ncchar('water');
name(nc{'mask_u'}.option_0_, 'option(0)')
name(nc{'mask_u'}.option_1_, 'option(1)')
%		 nc{'mask_u'}.FillValue_ = ncdouble(1);
 
nc{'mask_v'} = ncdouble('eta_v', 'xi_v'); %% 16770 elements.
nc{'mask_v'}.long_name = ncchar('mask on V-points');
nc{'mask_v'}.option_0_ = ncchar('land');
nc{'mask_v'}.option_1_ = ncchar('water');
name(nc{'mask_v'}.option_0_, 'option(0)')
name(nc{'mask_v'}.option_1_, 'option(1)')
%		 nc{'mask_v'}.FillValue_ = ncdouble(1);
 
nc{'mask_psi'} = ncdouble('eta_psi', 'xi_psi'); %% 16641 elements.
nc{'mask_psi'}.long_name = ncchar('mask on PSI-points');
nc{'mask_psi'}.option_0_ = ncchar('land');
nc{'mask_psi'}.option_1_ = ncchar('water');
name(nc{'mask_psi'}.option_0_, 'option(0)')
name(nc{'mask_psi'}.option_1_, 'option(1)')
%		 nc{'mask_psi'}.FillValue_ = ncdouble(1);

% Now, what about depths: "h" and "hraw".  <== DEPTHS.

% The following seems mistaken:
% [I think I fixed it: J Wilkin 2001/9/27]

nc{'angle'} = ncdouble('eta_rho', 'xi_rho'); %% 16900 elements.
% nc{'angle'}.long_name = ncchar('angle between xi axis and east');
nc{'angle'}.long_name = ncchar(...
    'angle ROMS xi axis is rotated anticlockwise from due east');
nc{'angle'}.equation = ncchar(...
    [ 'u(roms)=Real((Ueast+i*Vnorth)*exp(-i*angle)) and ' ...
      'v(roms)=Imag((Ueast+i*Vnorth)*exp(-i*angle))']); % JLW - for clarity
% nc{'angle'}.units = ncchar('degree');
nc{'angle'}.units = ncchar('radian'); % JLW


% Fill the variables with data.

disp(' ## Filling Variables...')

switch lower(projection)
case 'mercator'
		 theProjection = 'ME';
case 'stereographic'
		 theProjection = 'ST';
case 'lambert conformal conic'
		 theProjection = 'LC';
end

% Fill the variables.

% Need (x..., y...) in meters.  Currently, they are
%  in Mercator (projected) units.

nc{'JPRJ'}(:) = theProjection;
% nc{'spherical'}(:) = 'T';   % T or F -- uppercase okay?
% moved to bottom after file is closed
% nc{'xl'}(:) = xl;
% nc{'el'}(:) = el;

f = 2.*7.29e-5.*sin(geogrid_lat(j_rho, i_rho).*pi./180);
nc{'f'}(:) = f;

nc{'x_rho'}(:) = grid_x(j_rho, i_rho);
nc{'y_rho'}(:) = grid_y(j_rho, i_rho);

nc{'x_psi'}(:) = grid_x(j_psi, i_psi);
nc{'y_psi'}(:) = grid_y(j_psi, i_psi);

nc{'x_u'}(:) = grid_x(j_u, i_u);
nc{'y_u'}(:) = grid_y(j_u, i_u);

nc{'x_v'}(:) = grid_x(j_v, i_v);
nc{'y_v'}(:) = grid_y(j_v, i_v);

nc{'lon_rho'}(:) = geogrid_lon(j_rho, i_rho);
nc{'lat_rho'}(:) = geogrid_lat(j_rho, i_rho);

nc{'lon_psi'}(:) = geogrid_lon(j_psi, i_psi);
nc{'lat_psi'}(:) = geogrid_lat(j_psi, i_psi);

nc{'lon_u'}(:) = geogrid_lon(j_u, i_u);
nc{'lat_u'}(:) = geogrid_lat(j_u, i_u);

nc{'lon_v'}(:) = geogrid_lon(j_v, i_v);
nc{'lat_v'}(:) = geogrid_lat(j_v, i_v);

if ~isempty(bathymetry)
		 nc{'h'}(:) = bathymetry;
end

% Masking.

mask = ~~mask;
land = mask;
water = ~land;

rmask = water;

% Calculate other masking arrays.

umask = zeros(size(rmask));
vmask = zeros(size(rmask));
pmask = zeros(size(rmask));

for i = 2:LP
		 for j = 1:MP
		 		 umask(j, i-1) = rmask(j, i) * rmask(j, i-1);
		 end
end

for i = 1:LP
		 for j = 2:MP
		 		 vmask(j-1, i) = rmask(j, i) * rmask(j-1, i);
		 end
end

for i = 2:LP
		 for j = 2:MP
		 		 pmask(j-1, i-1) = rmask(j, i) * rmask(j, i-1) * rmask(j-1, i) * rmask(j-1, i-1);
		 end
end

nc{'mask_rho'}(:) = rmask;
nc{'mask_psi'}(:) = pmask(1:end-1, 1:end-1);
nc{'mask_u'}(:) = umask(1:end, 1:end-1);
nc{'mask_v'}(:) = vmask(1:end-1, 1:end);

% Average angle -- We should do this via (x, y) components.

temp = 0.5*(ang(1:end-1, :) + ang(2:end, :));
ang = zeros(n, m);
ang(2:2:end, 2:2:end) = temp;

nc{'angle'}(:) = ang(j_rho, i_rho);

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

% raw_grid_size = [m, n]

% geometry_sizes = [size(gx) size(gy)]

sx = 0.5*(gx(1:end-1, :) + gx(2:end, :));
sy = 0.5*(gy(:, 1:end-1) + gy(:, 2:end));

% raw_s_sizes = [size(sx) size(sy)]

% sx = sx(2:end-1, :);
% sy = sy(:, 2:end-1);

pm = 1 ./ sx;
pn = 1 ./ sy;

nc{'pm'}(:) = pm;
nc{'pn'}(:) = pn;

dmde = zeros(size(pm));
dndx = zeros(size(pn));

dmde(2:end-1, :) = 0.5*(1./pm(3:end, :) - 1./pm(1:end-2, :));
dndx(:, 2:end-1) = 0.5*(1./pn(:, 3:end) - 1./pn(:, 1:end-2));

nc{'dmde'}(:) = dmde;
nc{'dndx'}(:) = dndx;

% Final size of file:

s = size(nc);
disp([' ## Dimensions: ' int2str(s(1))])
disp([' ## Variables: ' int2str(s(2))])
disp([' ## Global Attributes: ' int2str(s(3))])
disp([' ## Record Dimension: ' name(recdim(nc))])
 
endef(nc)
close(nc)

nc_varput(theRomsFile,'spherical','T')
nc_varput(theRomsFile,'xl',xl)
nc_varput(theRomsFile,'el',el)
