function theResult = getbathymetry(self, theBathymetryFile)

% seagrid/getbathymetry -- Load and plot a bathymetry file.
%  getbathymetry(self, 'theBathymetryFile') loads and plots
%   the given bathymetry file on behalf of self, a "seagrid"
%   object.  If no filename is given, the current bathymetry
%   filename in self is used.  If a Mat-file, the variables
%   are expected to be "xbathy" (latitude), "ybathy" (longitude),
%   "zbathy" (arbitrary units, positive downwards).  If an
%   ascii file with three columns, the arrangement is expected
%   to be [xbathy ybathy zbathy].
 
% svn $Id$
%=======================================================================
% Copyright (C) 1999 Dr. Charles R. Denham, ZYDECO.
%  All Rights Reserved.
%   Disclosure without explicit written consent from the
%    copyright owner does not constitute publication.
%=======================================================================
 
% Version of 27-Apr-1999 08:48:25.
% Updated    11-Dec-2000 15:11:15.

if nargout > 0, theResult = self; end

if nargin < 1, help(mfilename), return, end

if nargin < 2
	theBathymetryFile = psget(self, 'itsBathymetryFile');
end

if ~isempty(theBathymetryFile) & ~any(theBathymetryFile == filesep)
	theBathymetryFile = which(theBathymetryFile);
end
psset(self, 'itsBathymetryFile', theBathymetryFile)

if isempty(theBathymetryFile), return, end

okay = 0;

xbathy = [];
ybathy = [];
zbathy = [];

s = load(theBathymetryFile);

switch class(s)
case 'double'
	if size(s, 2) == 3   % Three columns.
		xbathy = s(:, 1);
		ybathy = s(:, 2);
		zbathy = s(:, 3);
		okay = 1;
	end
case 'struct'
	try
		xbathy = s.xbathy;
		ybathy = s.ybathy;
		zbathy = s.zbathy;
		okay = 1;
	catch
	end
end

if ~okay
	disp([' ## Not a valid bathymetry file: ' theBathymetryFile])
	return
end

theBathymetryColor = psget(self, 'itsBathymetryColor');
theButtonDownFcn = get(gca, 'ButtonDownFcn');
h = findobj('Type', 'line', 'Tag', 'bathymetry');
if any(h), delete(h), end
theProjection = psget(self, 'itsProjection');
switch theProjection
case {'none', 'Geographic'}
	theProjection = 'Geographic';
	x = xbathy; y = ybathy;
otherwise
	sg_proj(theProjection)
	[x, y] = sg_ll2xy(xbathy, ybathy);
end
hold on
z = -zbathy;
h = plot3(x, y, z, '.', 'Color', theBathymetryColor, ...
					'MarkerSize', 10, 'Tag', 'bathymetry');
hold off
view(2)
zoomsafe 0
set([gca h], 'ButtonDownFcn', theButtonDownFcn)

self = doupdate(self);

if nargout > 0, theResult = self; end
