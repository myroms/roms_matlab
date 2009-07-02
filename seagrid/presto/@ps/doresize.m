function theResult = doresize(self, varargin)

% ps/doresize -- Process "resize" event for "ps".
%  doresize(self) processes a "resize" event on behalf
%   of self, a "ps" object.  All controls whose "UserData"
%   appears to be a "layout" are repositioned accordingly.

% svn $Id$
%=======================================================================
% Copyright (C) 1999 Dr. Charles R. Denham, ZYDECO.
%  All Rights Reserved.
%   Disclosure without explicit written consent from the
%    copyright owner does not constitute publication.
%=======================================================================
 
% Version of 02-Nov-1999 23:18:31.
% Updated    03-Nov-1999 01:29:18.

if nargout > 0, theResult = []; end
if nargin < 1, help(mfilename), return, end

theFigure = gcbf;

if isempty(theFigure)
	theFigure = handle(self);
	while ~isequal(get(theFigure, 'Type'), 'figure')
		theFigure = get(theFigure, 'Parent');
	end
end

oldUnits = get(theFigure, 'Units');
set(theFigure, 'Units', 'pixels')
theFigurePos = get(theFigure, 'Position');
set(theFigure, 'Units', oldUnits)

h = findobj(theFigure, 'Type', 'uicontrol');

for i = 1:length(h)
	u = get(h(i), 'UserData');
	if isstruct(u) & isfield(u, 'itsLayout')
		theLayout = u.itsLayout;
		if isequal(size(theLayout), [2 4])
			theNormalizedPos = theLayout(1, :);
			thePixelOffset = theLayout(2, :);
			thePos = theFigurePos([3:4 3:4]) .* ...
						theNormalizedPos + thePixelOffset;
			oldUnits = get(h(i), 'Units');
			set(h(i), 'Units', 'pixels')
			set(h(i), 'Position', thePos)
			set(h(i), 'Units', oldUnits)
		end
	end
end

if nargout > 0
	theResult = self;
end
