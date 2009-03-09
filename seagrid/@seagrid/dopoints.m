function theResult = dopoints(self)

% seagrid/dopoints -- Compute and draw points.
%  dopoints(self) draws the points associated
%   with self, a "seagrid" object.
 
% svn $Id$
%=======================================================================
% Copyright (C) 1999 Dr. Charles R. Denham, ZYDECO.
%  All Rights Reserved.
%   Disclosure without explicit written consent from the
%    copyright owner does not constitute publication.
%=======================================================================
 
% Version of 21-Apr-1999 09:03:33.
% Updated    28-Aug-2000 15:00:03.

if nargout > 0, theResult = self; end
if nargin < 1, help(mfilename), return, end

% Initialize.

thePoints = psget(self, 'itsPoints');

if length(thePoints)  < 4, return, end

x = zeros(size(thePoints));
y = zeros(size(thePoints));
theTags = cell(size(thePoints));
theCornerTag = psget(self, 'itsCornerTag');
theEdgePointTag = psget(self, 'itsEdgePointTag');
theCorners = [];   % Corner indices.
theEraseMode = psget(self, 'itsEraseMode');
theButtonDownFcn = 'psevent down';

thePointFlag = psget(self, 'itsPointFlag');

theCornerColor = psget(self, 'itsCornerColor');
theCornerMarker = psget(self, 'itsCornerMarker');
theEdgeColor = psget(self, 'itsEdgeColor');
theEdgeMarker = psget(self, 'itsEdgeMarker');

theMarkerSize = get(0, 'DefaultLineMarkerSize');
theMarkerSize = 15;
theMarkerSize = 9;

% Get the xy data and corner-point indices.

[m, n] = size(thePoints);
if min(size(thePoints)) == 1
	for k = 1:length(thePoints)
		x(k) = get(thePoints(k), 'XData');
		y(k) = get(thePoints(k), 'YData');
		theTags{k} = get(thePoints(k), 'Tag');
		switch theTags{k}
		case 'corner-point'
			theCorners = [theCorners k];   % #1 is always a theCorners-point.
		otherwise
		end
	end
else
	x = thePoints(:, 1); x = x(:).';
	y = thePoints(:, 2); y = y(:).';
	t = thePoints(:, 3); t = t(:).';
	theCorners = find(t);
	theTags = cell(size(t));
	for k = 1:length(t)
		switch t(k)
		case 0
			theTags{k} = theEdgePointTag;
		otherwise
			theTags{k} = theCornerTag;
		end
	end
	thePoints = zeros(size(x));
end

delete(findobj('Type', 'line', 'Tag', theEdgePointTag))
delete(findobj('Type', 'line', 'Tag', theCornerTag))

% Delete old points.  We need to recreate the points
%  to keep the screen from flashing, which would occur
%  if we were to use "bringtofront" instead.

% Draw new points.

thePoints = zeros(1, length(x));
for k = 1:length(x)
	hold on
	thePoints(k) = plot(x(k), y(k), 'EraseMode', theEraseMode);
	hold off
	theColor = theEdgeColor;
	theMarker = theEdgeMarker;
	if any(k == theCorners)
		theColor = theCornerColor;
		theMarker = theCornerMarker;
		if k == 1, theMarker = '*'; end
	end
	set(thePoints(k), 'Marker', theMarker, ...
			'Color', theColor, ...
			'Tag', theTags{k})
end
	
set(thePoints(theCorners), 'MarkerFaceColor', theColor)
		
psset(self, 'itsPoints', thePoints)

switch thePointFlag
case 0
	theVisible = 'off';
otherwise
	theVisible = 'on';
end

set(thePoints, 'ButtonDownFcn', theButtonDownFcn, ...
				'LineStyle', 'none', 'MarkerSize', theMarkerSize, ...
				'UserData', thePoints, 'Visible', theVisible)

if nargout > 0, theResult = self; end
