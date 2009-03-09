function theResult = dosave(self, theFilename)

% seagrid/dosave -- Save a grid.
%  dosave(self, 'theFilename') saves the grid and other
%   information associated with self, a "seagrid" object,
%   to 'theFilename'.  The "uiputfile" dialog is invoked
%   if an empty or wildcarded name is provided.  A grid
%   if size M-by-N contains (M+1)-by-(N+1) points that
%   describe the cell corners.  The bathymetric and
%   masking points lie at the cell-centers, so their
%   dimension is M-by-N.
 
% svn $Id$
%=======================================================================
% Copyright (C) 1999 Dr. Charles R. Denham, ZYDECO.
%  All Rights Reserved.
%   Disclosure without explicit written consent from the
%    copyright owner does not constitute publication.
%=======================================================================
 
% Version of 13-May-1999 15:58:14.
% Updated    11-Dec-2000 15:16:12.

RCF = 180 / pi;

if nargout > 0, theResult = self; end
if nargin < 1, help(mfilename), return, end
if nargin < 2, theFilename = ''; end

if nargout > 0, theResult = self; end

if isempty(theFilename), theFilename = '*.mat'; end

if any(theFilename == '*')
	theSuggested = 'seagrid.mat';
	[theFile, thePath] = uiputfile(theSuggested, 'Save As SeaGrid File:');
	if ~any(theFile), return, end
	if thePath(end) ~= filesep, thePath(end+1) = filesep; end
	theFilename = [thePath theFile];
end

psset(self, 'itsSeaGridOutputFile', theFilename)

% Double the grid-size -- no longer.
%  We will let "seagrid2ecom" and other
%  utilities take care of this by simple
%  interpolation.

if (0)
	theOldGridSize = psget(self, 'itsGridSize');
	theTempGridSize = 2*theOldGridSize;
	psset(self, 'itsGridSize', theTempGridSize);
	if (0)
		dospacings(self)
		dogrid(self, 1)
	else
		doupdate(self, 1)
	end
end

% Turn masktool off.

theMaskToolFlag = psget(self, 'itsMaskToolFlag');
theMaskTool = psget(self, 'itsMaskTool');
if any(theMaskTool)
	domasktool(self)   % Turn it off.
end

s.created_on = datestr(now);
s.created_by = mfilename;

s.projection = psget(self, 'itsProjection');
s.projection_center = psget(self, 'itsProjectionCenter');
s.longitude_bounds = psget(self, 'itsLongitudeBounds');
s.latitude_bounds = psget(self, 'itsLatitudeBounds');

s.bathymetry = psget(self, 'itsBathymetryFile');
s.coastline = psget(self, 'itsCoastlineFile');

theCornerTag = psget(self, 'itsCornerTag');

h = psget(self, 'itsPoints');
x = zeros(size(h));
y = zeros(size(h));
t = zeros(size(h));

for k = 1:length(h)
	x(k) = get(h(k), 'XData');
	y(k) = get(h(k), 'YData');
	if isequal(get(h(k), 'Tag'), theCornerTag)
		t(k) = 1;
	end
end
s.points = [x(:) y(:) t(:)];

s.end_slope_flag = psget(self, 'itsEndSlopeFlag');

theGridSize = psget(self, 'itsGridSize');

s.grid_size = theGridSize;

s.grids = psget(self, 'itsGrids');

s.spacings = psget(self, 'itsSpacings');
s.default_spacings = psget(self, 'itsDefaultSpacings');
s.spaced_edges = psget(self, 'itsSpacedEdges');
s.spaced_grids = psget(self, 'itsSpacedGrids');

u = s.spaced_grids{1}; v = s.spaced_grids{2};
sg_proj(s.projection)
[lon, lat] = sg_xy2ll(u, v);
s.geographic_grids = {lon, lat};

b = -psget(self, 'itsGriddedBathymetry');   % Note negative.
c = psget(self, 'itsClippingDepths');

s.gridded_bathymetry = b;
s.clipping_depths = c;

if (0), b = b(2:2:end-1, 2:2:end-1); end   % No longer.

b(b > c(2)) = c(2);
b(b < c(1)) = c(1);

s.bottom = b;
s.top = zeros(size(b)) + c(1);

dx = earthdist(lon(:, 2:end), lat(:, 2:end), lon(:, 1:end-1), lat(:, 1:end-1));
dy = earthdist(lon(2:end, :), lat(2:end, :), lon(1:end-1, :), lat(1:end-1, :));

s.geometry = {dx, dy};

% Grid-cell orientation, degrees counter-clockwise from
%  east, presently based on flat-earth approximation.

dlon = diff(lon.').';
dlat = diff(lat.').';
clat = cos(lat / RCF);
clat(:, end) = [];
s.orientation = atan2(dlat, dlon .* clat) * RCF;

theMask = psget(self, 'itsMask');
theWater = psget(self, 'itsWater');
theLand = psget(self, 'itsLand');
% theMask = theMask(1:end-1, 1:end-1)

s.mask = theMask;
s.water = theWater;
s.land = theLand;

save(theFilename, 's')

% Turn masktool on.

if any(theMaskToolFlag)
	f = psget(self, 'itsPoints');
	if any(f), set(f, 'Visible', 'off'), end
	f = psget(self, 'itsEdges');
	if any(f), set(f, 'Visible', 'off'), end
	domasktool(self);   % Turn it on.
end

% Revert to original grid-size -- not needed anymore.

if (0)
	psset(self, 'itsGridSize', theOldGridSize);
	if (0)
		dospacings(self)
		dogrid(self, 1)
	else
		doupdate(self, 1)
	end
end

if nargout > 0, theResult = self; end
