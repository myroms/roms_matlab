function seagrid2scrum(varargin)

% function seagrid2scrum(theSeagridFile, theScrumFile, theGridTitle)

% seagrid2scrum -- Output SCRUM format from SeaGrid file.
%  seagrid2scrum('theSeagridFile', 'theScrumFile', 'theGridTitle')
%   creates a SCRUM file, based on the given SeaGrid file.  The
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
% Updated    31-May-2000 10:33:39.

if (1)
	disp([' ## seagrid2scrum is obsolete; calling "seagrid2roms" instead.'])
	if nargin < 1
		feval('seagrid2roms')
	else
		feval('seagrid2roms', varargin{:})
	end
	return
end

if nargin < 1, theSeagridFile = '*.mat'; end
if nargin < 2, theScrumFile = 'scrum_model_grid.nc'; end
if nargin < 3, theGridTitle = char(zeros(1, 128)+abs(' ')); end

if isempty(theSeagridFile) | any(theSeagridFile == '*')
	[f, p] = uigetfile(theSeagridFile, 'Select SeaGrid File:');
	if ~any(f), return, end
	if p(end) ~= filesep, p(end+1) = filesep; end
	theSeagridFile = [p f]
end

if nargin < 2 | isempty(theSeagridFile) | any(theScrumFile == '*')
	[f, p] = uiputfile(theScrumFile, 'Save to Scrum File:');
	if ~any(f), return, end
	if p(end) ~= filesep, p(end+1) = filesep; end
	theScrumFile = [p f]
end

disp([' ## SeaGrid Source File  : ' theSeagridFile])
disp([' ## Scrum Destination File: ' theScrumFile])

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

grid_x = s.grids{1};
grid_y = s.grids{2};
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

bathymetry = s.gridded_bathymetry;
projection = s.projection;
ang = s.orientation;
min_depth = s.clipping_depths(1);
max_depth = s.clipping_depths(2);
coriolis = geogrid_lon;

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

[m, n] = size(grid_x);

% Flip arrays top for bottom.
%  Is this necessary?

FLIPPING = 0;
FLIPPING = 1;

if FLIPPING
	grid_x = flipud(grid_x);
	grid_y = flipud(grid_y);
	geogrid_lon = flipud(geogrid_lon);
	geogrid_lat = flipud(geogrid_lat);
	geometry{1} = flipud(geometry{1});
	mask = flipud(mask);
	bathymetry = flipud(bathymetry);
	ang = flipud(ang);
	coriolis = flipud(coriolis);
	spaced_x = flipud(spaced_x);
	spaced_y = flipud(spaced_y);
end

% Create the Scrum NetCDF file.

%% ncdump('priapus:MATLAB:toolbox:scrum:tasman_grd.nc')   %% Generated 17-Jun-1999 17:37:30
 
nc = netcdf(theScrumFile, 'clobber');
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
%  sub-grids for the SCRUM scheme, ignoring the
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
i_u = 2:2:m-1; j_u = 3:2:n-2;
i_v = 3:2:m-2; j_v = 2:2:n-1;

% The xi direction (left-right):
%	LP = (m-1)/2;   % The rho dimension.
LP = (n-1)/2;   % The rho dimension.
L = LP-1;       % The psi dimension.

% The eta direction (up-down):
%	MP = (n-1)/2;   % The rho dimension.
MP = (m-1)/2;   % The rho dimension.
M = MP-1;       % The psi dimension.

[L M LP MP]

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
 
nc{'f0'} = ncdouble; %% 1 element.
nc{'f0'}.long_name = ncchar('Coriolis parameter central value on a beta-plane');
%	nc{'f0'}.FillValue_ = ncdouble(0);
 
nc{'dfdy'} = ncdouble; %% 1 element.
nc{'dfdy'}.long_name = ncchar('Coriolis parameter gradient on a beta-plane');
%	nc{'dfdy'}.FillValue_ = ncdouble(0);
 
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
%	nc{'mask_u'}.FillValue_ = ncdouble(1);
 
nc{'mask_v'} = ncdouble('eta_v', 'xi_v'); %% 16770 elements.
nc{'mask_v'}.long_name = ncchar('mask on V-points');
nc{'mask_v'}.option_0_ = ncchar('land');
nc{'mask_v'}.option_1_ = ncchar('water');
name(nc{'mask_v'}.option_0_, 'option(0)')
name(nc{'mask_v'}.option_1_, 'option(1)')
%	nc{'mask_v'}.FillValue_ = ncdouble(1);
 
nc{'mask_psi'} = ncdouble('eta_psi', 'xi_psi'); %% 16641 elements.
nc{'mask_psi'}.long_name = ncchar('mask on PSI-points');
nc{'mask_psi'}.option_0_ = ncchar('land');
nc{'mask_psi'}.option_1_ = ncchar('water');
name(nc{'mask_psi'}.option_0_, 'option(0)')
name(nc{'mask_psi'}.option_1_, 'option(1)')
%	nc{'mask_psi'}.FillValue_ = ncdouble(1);

% Now, what about depths: "h" and "hraw".  <== DEPTHS.

% The following seems mistaken:
 
nc{'angle'} = ncdouble('eta_rho', 'xi_rho'); %% 16900 elements.
nc{'angle'}.long_name = ncchar('angle between xi axis and east');
nc{'angle'}.units = ncchar('degree');

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

nc{'JPRJ'}(:) = theProjection;
nc{'spherical'}(:) = 'T';

nc{'x_rho'}(:) = grid_x(i_rho, j_rho);
nc{'y_rho'}(:) = grid_y(i_rho, j_rho);

nc{'x_psi'}(:) = grid_x(i_psi, j_psi);
nc{'y_psi'}(:) = grid_y(i_psi, j_psi);

% [size(nc{'x_u'}) size(grid_x(i_u, j_u))]

nc{'x_u'}(:) = grid_x(i_u, j_u);
nc{'y_u'}(:) = grid_y(i_u, j_u);

nc{'x_v'}(:) = grid_x(i_v, j_v);
nc{'y_v'}(:) = grid_y(i_v, j_v);

nc{'lon_rho'}(:) = geogrid_lon(i_rho, j_rho);
nc{'lat_rho'}(:) = geogrid_lat(i_rho, j_rho);

nc{'lon_psi'}(:) = geogrid_lon(i_psi, j_psi);
nc{'lat_psi'}(:) = geogrid_lat(i_psi, j_psi);

nc{'lon_u'}(:) = geogrid_lon(i_u, j_u);
nc{'lat_u'}(:) = geogrid_lat(i_u, j_u);

nc{'lon_v'}(:) = geogrid_lon(i_v, j_v);
nc{'lat_v'}(:) = geogrid_lat(i_v, j_v);

whos

nc{'h'}(:) = bathymetry;

water = zeros(m, n);
water(2:2:end, 2:2:end) = 1 - mask;

nc{'mask_rho'}(:) = water(i_rho, j_rho);
nc{'mask_psi'}(:) = water(i_psi, j_psi);
nc{'mask_u'}(:) = water(i_u, j_u);
nc{'mask_v'}(:) = water(i_v, j_v);

temp = 0.5*(ang(1:end-1, :) + ang(2:end, :));
ang = zeros(m, n);
ang(2:2:end, 2:2:end) = temp;

nc{'angle'}(:) = ang(i_rho, j_rho);

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

% p_sizes = [size(pm) size(pn)]

dmde = zeros(size(pm));
dndx = zeros(size(pn));

dmde(:, 2:end-1) = 0.5*(1./pm(:, 3:end) - 1./pm(:, 1:end-2));
dndx(2:end-1, :) = 0.5*(1./pn(3:end, :) - 1./pn(1:end-2, :));

% d_sizes = [size(dmde) size(dndx)]

% dmde = dmde(:, 1:2:end);
% dndx = dndx(1:2:end, :);

nc{'pm'}(:) = pm;
nc{'pn'}(:) = pn;

% [size(dmde) size(nc{'dmde'})]

nc{'dmde'}(:) = dmde;

% [size(dndx) size(nc{'dndx'})]

nc{'dndx'}(:) = dndx;

% Final size of file:

s = size(nc);
disp([' ## Dimensions: ' int2str(s(1))])
disp([' ## Variables: ' int2str(s(2))])
disp([' ## Global Attributes: ' int2str(s(3))])
disp([' ## Record Dimension: ' name(recdim(nc))])
 
endef(nc)
close(nc)
