function theResult = doload(self, theFilename)

% seagrid/doload -- Load a grid.
%  doload(self, 'theFilename') loads the grid and other
%   information associated with self, a "seagrid" object,
%   from 'theFilename'.  The "uigetfile" dialog is invoked
%   if an empty or wildcarded name is provided.
 
% svn $Id$
%=======================================================================
% Copyright (C) 1999 Dr. Charles R. Denham, ZYDECO.
%  All Rights Reserved.
%   Disclosure without explicit written consent from the
%    copyright owner does not constitute publication.
%=======================================================================
 
% Version of 13-May-1999 15:58:14.
% Updated    10-Jul-2000 10:07:06.

if nargout > 0, theResult = []; end
if nargin < 1, help(mfilename), return, end
if nargin < 2, theFilename = ''; end

if nargout > 0, theResult = self; end

if isempty(theFilename), theFilename = '*.mat'; end

if any(theFilename == '*')
	[theFile, thePath] = uigetfile(theFilename, 'Select SeaGrid File:');
	if ~any(theFile), return, end
	if thePath(end) ~= filesep, thePath(end+1) = filesep; end
	theFilename = [thePath theFile];
end

loaded = load(theFilename);
if ~isfield(loaded, 's'), return, end

theMaskToolFlag = psget(self, 'itsMaskToolFlag');
theMaskTool = psget(self, 'itsMaskTool');
if any(theMaskTool)
	domasktool(self)   % Turn it off.
end

s = loaded.s;

psset(self, 'itsSeaGridOutputFile', theFilename)
psset(self, 'itsProjection', s.projection)
psset(self, 'itsProjectionCenter', s.projection_center)
psset(self, 'itsLongitudeBounds', s.longitude_bounds)
psset(self, 'itsLatitudeBounds', s.latitude_bounds)
psset(self, 'itsCoastlineFile', s.coastline)
psset(self, 'itsBathymetryFile', s.bathymetry)
psset(self, 'itsPoints', s.points)
psset(self, 'itsEndSlopeFlag', s.end_slope_flag)
psset(self, 'itsGridSize', s.grid_size)
psset(self, 'itsGrids', s.grids)
psset(self, 'itsSpacings', s.spacings)
psset(self, 'itsDefaultSpacings', s.default_spacings)
psset(self, 'itsSpacedEdges', s.spaced_edges)
psset(self, 'itsSpacedGrids', s.spaced_grids)
psset(self, 'itsGriddedBathymetry', -s.gridded_bathymetry)
psset(self, 'itsClippingDepths', s.clipping_depths)
psset(self, 'itsMask', s.mask)
psset(self, 'itsWater', s.water)
psset(self, 'itsLand', s.land)

sg_proj(s.projection)

if (0), psset(self, 'itsGridSize', s.grid_size/2), end   % No longer.

getboundary(self, 'clear')

% dogrid(self, 1)

% If alive, the masktool will be turned back on during update.

doupdate(self, 1)

getcoastline(self)
getbathymetry(self)

if nargout > 0, theResult = self; end
