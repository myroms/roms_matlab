function theResult = dosetup(self)

% seagrid/dosetup -- Dialog for "seagrid" parameters.
%  setup(self) presents a dialog for self, a "seagrid"
%   object.
 
% svn $Id$
%=======================================================================
% Copyright (C) 1999 Dr. Charles R. Denham, ZYDECO.
%  All Rights Reserved.
%   Disclosure without explicit written consent from the
%    copyright owner does not constitute publication.
%=======================================================================
 
% Version of 07-Apr-1999 22:15:27.
% Updated    10-Jan-2000 16:30:01.

if nargin < 1, help(mfilename), return, end

busy

if (0)
	theProjections = ...
		{'Geographic', 'Lambert Conformal Conic', 'Mercator', 'Stereographic'};
end

theProjections = ...
	{'Mercator', 'Stereographic', '-', 'Lambert Equal Area'};
	
theProjection = psget(self, 'itsProjection');

proj_code = 1;
for i = 1:length(theProjections)
	if isequal(lower(theProjection), lower(theProjections{i}))
		proj_code = i;
		break;
	end
end

theProjectionCenter = psget(self, 'itsProjectionCenter');
if isempty(theProjectionCenter)
	theProjectionCenter = [0 0 0];
	psset(self, 'itsProjectionCenter', theProjectionCenter)
end
theLongitudeBounds = psget(self, 'itsLongitudeBounds');
theLatitudeBounds = psget(self, 'itsLatitudeBounds');

theGridSize = psget(self, 'itsGridSize');
if length(theGridSize) == 1
	theGridSize = theGridSize .* [1 1];
end

theGridDensityFactor = psget(self, 'itsGridDensityFactor');

theClippingDepths = psget(self, 'itsClippingDepths');

theEndSlopeFlag = psget(self, 'itsEndSlopeFlag');

useMexFile = psget(self, 'itUsesMexFile');

s = [];

s.Projection = {theProjections, proj_code};
s.Proj.CentralLongitude = theProjectionCenter(1);
s.Proj.CentralLatitude = theProjectionCenter(2);
s.Proj.RotationAngle = theProjectionCenter(3);
s.Proj.LongitudeMinMax = theLongitudeBounds;
s.Proj.LatitudeMinMax = theLatitudeBounds;

theGridSize = max(theGridSize, [6 6]);   % Four cells minimum.
theCellCount = theGridSize;
s.Cells_Edge_1 = theCellCount(1);
s.Cells_Edge_2 = theCellCount(2);
s.MinimumDepth = theClippingDepths(1);
s.MaximumDepth = theClippingDepths(2);

s.EndSlopeFlag = {'checkbox', any(theEndSlopeFlag)};

s.UseMexFile = {'checkbox', any(useMexFile)};

SeaGrid_Setup = s;

if (0)
%	s = uigetinfo(SeaGrid_Setup);
else
	s = guido(SeaGrid_Setup);
end

if ~isempty(s)
	theNewProjection = getinfo(s, 'Projection');
	if isequal(theNewProjection, '-')
		theNewProjection = theProjection;
	end
	t = s.Proj;
	theNewCenterLon = getinfo(t, 'CentralLongitude');
	theNewCenterLat = getinfo(t, 'CentralLatitude');
	theNewCenterAng = getinfo(t, 'RotationAngle');
	theNewProjectionCenter = [theNewCenterLon theNewCenterLat theNewCenterAng];
	psset(self, 'itsNewProjection', theNewProjection)
	psset(self, 'itsNewProjectionCenter', theNewProjectionCenter)
	
	doreproject(self)
	
	theLongitudeBounds = getinfo(s, 'Proj.LongitudeMinMax');
	theLatitudeBounds = getinfo(s, 'Proj.LatitudeMinMax');
	theCells_Edge_1 = getinfo(s, 'Cells_Edge_1');
	theCells_Edge_2 = getinfo(s, 'Cells_Edge_2');
	theGridLines_Edge_1 = theCells_Edge_1;
	theGridLines_Edge_2 = theCells_Edge_2;
	psset(self, 'itsLongitudeBounds', theLongitudeBounds)
	psset(self, 'itsLatitudeBounds', theLatitudeBounds)
	if ~isempty(theGridLines_Edge_1) & ~isempty(theGridLines_Edge_2)
		theGridSize = [theGridLines_Edge_1 theGridLines_Edge_2];
		theGridSize = ceil(max(theGridSize, [6 6]));
		psset(self, 'itsGridSize', theGridSize)
	end
	theGridDensityFactor = getinfo(s, 'GridDensityFactor');
	if 0 & ~isempty(theGridDensityFactor)
		psset(self, 'itsGridDensityFactor', theGridDensityFactor)
	end
	theClippingDepths = getinfo(s, 'ClippingDepths');
	if 0 & ~isempty(theClippingDepths)
		psset(self, 'itsClippingDepths', theClippingDepths);
	end
	theMinimumDepth = getinfo(s, 'MinimumDepth');
	theMaximumDepth = getinfo(s, 'MaximumDepth');
	if ~isempty(theMinimumDepth) & ~isempty(theMaximumDepth)
		theClippingDepths = [theMinimumDepth theMaximumDepth];
		psset(self, 'itsClippingDepths', theClippingDepths);
	end
	theEndSlopeFlag = getinfo(s, 'EndSlopeFlag');
	if ~isempty(theEndSlopeFlag)
		psset(self, 'itsEndSlopeFlag', any(theEndSlopeFlag))
	end
	useMexFile = getinfo(s, 'UseMexFile');
	if ~isempty(useMexFile)
		psset(self, 'itUsesMexFile', any(useMexFile));
	end
	psset(self, 'itNeedsUpdate', 1)
else
	psset(self, 'itNeedsUpdate', 0)
end

if nargout > 0, theResult = self; end

idle
