function theResult = dosave(self, theFilename)

% seagrid/dosave -- Save a grid.
%  dosave(self, 'theFilename') saves the grid and other
%   information associated with self, a "seagrid" object,
%   to 'theFilename'.  The "uiputfile" dialog is invoked
%   if an empty or wildcarded name is provided.
 
% svn $Id$
%=======================================================================
% Copyright (C) 1999 Dr. Charles R. Denham, ZYDECO.
%  All Rights Reserved.
%   Disclosure without explicit written consent from the
%    copyright owner does not constitute publication.
%=======================================================================
 
% Version of 13-May-1999 15:58:14.
% Updated    28-Mar-2000 10:17:18.

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

% Double the grid-size.

if (1)
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

b = b(2:2:end-1, 2:2:end-1);
b(b > c(2)) = c(2);
b(b < c(1)) = c(1);

s.bottom = b;
s.top = zeros(size(b)) + c(1);

dx = earthdist(lon(:, 2:end), lat(:, 2:end), lon(:, 1:end-1), lat(:, 1:end-1));
dy = earthdist(lon(2:end, :), lat(2:end, :), lon(1:end-1, :), lat(1:end-1, :));

s.geometry = {dx, dy};

% Grid-cell orientation counter-clockwise from east,
%  presently based on flat-earth approximation.

dlon = diff(lon.').';
dlat = diff(lat.').';
clat = cos(lat / RCF);
clat(:, end) = [];
s.orientation = atan2(dlat, dlon .* clat) * RCF;

s.mask = psget(self, 'itsMask');

save(theFilename, 's')

% Revert to original grid-size.

if (1)
	psset(self, 'itsGridSize', theOldGridSize);
	if (0)
		dospacings(self)
		dogrid(self, 1)
	else
		doupdate(self, 1)
	end
end

if nargout > 0, theResult = self; end
