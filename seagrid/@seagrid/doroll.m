function theResult = doroll(self)

% seagrid/doroll -- Roll the main corner.
%  doroll(self) rolls the main grid corner counter-
%   clockwise, on behalf of self, a "seagrid" object.
 
% svn $Id$
%=======================================================================
% Copyright (C) 1999 Dr. Charles R. Denham, ZYDECO.
%  All Rights Reserved.
%   Disclosure without explicit written consent from the
%    copyright owner does not constitute publication.
%=======================================================================
 
% Version of 26-Apr-1999 10:29:08.
% Updated    21-Sep-2000 15:05:25.

if nargout > 0, theResult = self; end
if nargin < 1, help(mfilename), return, end

thePoints = psget(self, 'itsPoints');
theCorners = [];
theCornerTag = psget(self, 'itsCornerTag');
for i = 1:length(thePoints)
	theTag = get(thePoints(i), 'Tag');
	if isequal(theTag, theCornerTag)
		theCorners = [theCorners i];
	end
end
i = [theCorners(2):length(thePoints) theCorners(1):theCorners(2)-1];
thePoints = thePoints(i);
psset(self, 'itsPoints', thePoints)
theGridSize = psget(self, 'itsGridSize');
psset(self, 'itsGridSize', theGridSize([2 1]))

theSpacings = psget(self, 'itsSpacings');
u = theSpacings{1}; u = 1 - u(end:-1:1);
v = theSpacings{2};
psset(self, 'itsSpacings', {v u})
theGrids = psget(self, 'itsGrids');
u = theGrids{1}; u = flipud(fliplr(u.'));
v = theGrids{2}; v = flipud(fliplr(v.'));
psset(self, 'itsGrids', {u v})

doupdate(self, 1)

if nargout > 0, theResult = self; end
