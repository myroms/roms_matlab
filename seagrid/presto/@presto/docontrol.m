function theResult = docontrol(self, theEvent, varargin)

% presto/docontrol -- Handler for "presto" controls.
%  docontrol(self, 'theEvent') handles 'theEvent' on
%   behalf of self, a "presto" object.
 
% svn $Id$
%=======================================================================
% Copyright (C) 1999 Dr. Charles R. Denham, ZYDECO.
%  All Rights Reserved.
%   Disclosure without explicit written consent from the
%    copyright owner does not constitute publication.
%=======================================================================
 
% Version of 05-Nov-1999 11:49:51.
% Updated    05-Nov-1999 15:00:06.

if nargout > 0, theResult = []; end
if nargin < 1, help(mfilename), return, end

switch theEvent
case {'bottom', 'right', 'top', 'left'}
	theValue = get(gcbo, 'Value');
	theMin = get(gcbo, 'Min');
	theMax = get(gcbo, 'Max');
	disp([' ## ' theEvent ' ' mat2str([theMin theValue theMax])])
otherwise
	theValue = get(gcbo, 'Value');
	disp([' ## ' theEvent mat2str(theValue)])
end

if nargout > 0
	theResult = self;
end
