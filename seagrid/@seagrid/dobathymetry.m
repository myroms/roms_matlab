function theResult = dobathymetry(self, itNeedsUpdate)

% seagrid/dobathymetry -- Attach bathymetry to grid.
%  dobathymetry(self) grids and draws bathymetry on
%   the grid, on behalf of self, a "seagrid" object.
%   Only the grid-cell centers are represented.
 
% svn $Id$
%=======================================================================
% Copyright (C) 1999 Dr. Charles R. Denham, ZYDECO.
%  All Rights Reserved.
%   Disclosure without explicit written consent from the
%    copyright owner does not constitute publication.
%=======================================================================
 
% Version of 27-Apr-1999 16:40:37.
% Updated    15-Nov-2000 15:20:22.

if nargout > 0, theResult = self; end
if nargin < 1, help(mfilename), return, end
if nargin < 2, itNeedsUpdate = psget(self, 'itNeedsUpdate'); end

theFigure = ps(self);
isVerbose = psget(self, 'itIsVerbose');
theGridSize = psget(self, 'itsGridSize');

h = findobj(theFigure, 'Type', 'line', 'Tag', 'bathymetry');
if ~any(h)
	if nargout > 0, theResult = self; end
	return
end

delete(colorbar)

x = get(h, 'XData');
y = get(h, 'YData');
z = get(h, 'ZData');

if isempty(x) | isempty(y) | isempty(z)
	if nargout > 0, theResult = self; end
	return
end

% Bathymetry at grid-cell centers only.

theSpacedGrids = psget(self, 'itsSpacedGrids');
if isempty(theSpacedGrids), return, end

xi = theSpacedGrids{1};   % Corners.
oldxi = xi;
xi = interp2(xi, 1);
xi = xi(2:2:end-1, 2:2:end-1);   % Centers.

yi = theSpacedGrids{2};   % Corners.
oldyi = yi;
yi = interp2(yi, 1);
yi = yi(2:2:end-1, 2:2:end-1);   % Centers.

zi = psget(self, 'itsGriddedBathymetry');

if ~isequal(size(xi), size(yi)) | ~isequal(size(xi), size(zi))
	itNeedsUpdate = 1;
end

% We grid the bathymetry as best we can, then re-grid
%  the NaN results using the "nearest" depth data.

if itNeedsUpdate
	busy
	drawnow
	if isunix
		shouldAlert = (prod(theGridSize) > 1000);
	else
		shouldAlert = (prod(theGridSize) > 100);
	end
	if shouldAlert
		h = warndlg('Please wait ...', 'Computing Depths');
		drawnow
	end
	if (1)
		theGriddingMethod = psget(self, 'itsGriddingMethod');
		zi = griddata(x, y, z, xi, yi, theGriddingMethod);   % Grid z.
	else
		theGriddingMethod = [inf, 1, 1];
		zi = griddata1(x, y, z, xi, yi, theGriddingMethod);   % Grid z.
	end
	f = find(isnan(zi));
	
	if isVerbose
		disp([' ## ' mfilename ': NaNs not modified.'])
	end
	f = [];   % Testing only -- do not modify NaNs.

	if any(f)   % Use "nearest" depth point.
		tri = delaunay(x, y);
		indices = dsearch(x, y, tri, xi(f), yi(f));
		if all(indices)
			zi(f) = z(indices);
		end
	end
	psset(self, 'itsGriddedBathymetry', zi);
	if shouldAlert & ishandle(h), delete(h), end
	idle
else
	zi = psget(self, 'itsGriddedBathymetry');
end

% Mask the bathymetry.

theMask = psget(self, 'itsMask');
if ~isempty(theMask)
	if (0)   % Old way.
		b = zi(2:2:end, 2:2:end);
	else
		b = zi;
	end
	b(logical(theMask)) = NaN;
	if (0)   % old way.
		zi(2:2:end, 2:2:end) = b;
	else
		zi = b;
	end
	psset(self, 'itsGriddedBathymetry', zi);
end

theFigure = ps(self);
setsafe(0, 'CurrentFigure', theFigure)

hold on

% Note: we will want to expand the (xi, yi, zi) by
%  one row and column to fill the drawn grid, for
%  appearances' sake.

if (0)
	zi = zi;
else
	newzi = zi;
	newzi(end+1, :) = newzi(end, :);
	newzi(:, end+1) = newzi(:, end);
end

if (0)
	
if (0)
	h = surf(xi, yi, zi);
else
	h = surf(oldxi, oldyi, newzi);
end

theEraseMode = 'xor';

theHitTest = 'off';
theHitTest = 'on';

set(h, ...
	'EraseMode', theEraseMode, ...
	'HitTest', theHitTest, ...
	'Tag', 'gridded-bathymetry' ...
	);
	
if (0), colormap copper, end   % Color needs work
shading flat
view(2)

tic
if (0)
	[c, h] = contour(xi, yi, zi);
else
	[c, h] = contour(xi, yi, zi);
end
set(h, ...
		'EraseMode', theEraseMode, ...
		'HitTest', theHitTest, ...
		'Tag', 'contoured-bathymetry' ...
	);
t = toc;

hold off

h = colorbar;
setsafe(get(h, 'Ylabel'), 'String', 'Depth')


end

if nargout > 0, theResult = self; end
