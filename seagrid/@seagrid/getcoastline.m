function theResult = getcoastline(self, theCoastlineFile)

% seagrid/getcoastline -- Load and plot a coastline.
%  getcoastline(self) loads and plots a coastline file
%   on behalf of self, a "seagrid" object.  If no
%   filename is given, the current coastline filename
%   in self is used.
 
% svn $Id$
%=======================================================================
% Copyright (C) 1999 Dr. Charles R. Denham, ZYDECO.
%  All Rights Reserved.
%   Disclosure without explicit written consent from the
%    copyright owner does not constitute publication.
%=======================================================================
 
% Version of 27-Apr-1999 08:48:25.
% Updated    11-Dec-2000 15:11:40.

if nargout > 0, theResult = self; end

if nargin < 1, help(mfilename), return, end

if nargin < 2
	theCoastlineFile = psget(self, 'itsCoastlineFile');
end

if ~isempty(theCoastlineFile) & ~any(theCoastlineFile == filesep)
	theCoastlineFile = which(theCoastlineFile);
end
psset(self, 'itsCoastlineFile', theCoastlineFile)

if isempty(theCoastlineFile), return, end

load(theCoastlineFile)

theCoastlineColor = psget(self, 'itsCoastlineColor');
theButtonDownFcn = get(gca, 'ButtonDownFcn');
h = findobj('Type', 'line', 'Tag', 'coastline');
if any(h), delete(h), end

theProjection = psget(self, 'itsProjection');
switch theProjection   % Needs cleaning up.
case {'none', 'Geographic'}
	theProjection = 'Geographic';
	x = lon; y = lat;
otherwise
	sg_proj(theProjection)
	[x, y] = sg_ll2xy(lon, lat);
	x = real(x);
	y = real(y);
end

hold on
h = plot(x, y, '-', 'Color', theCoastlineColor, ...
						'Tag', 'coastline');
hold off
zoomsafe 0
set([gca h], 'ButtonDownFcn', theButtonDownFcn)
switch theProjection
case 'none'
	xlabel('Longitude')
	ylabel('Latitude')
otherwise
	xlabel('X')
	ylabel('Y')
	title([theProjection ' Projection'])
end

self = doupdate(self);

if nargout > 0, theResult = self; end
