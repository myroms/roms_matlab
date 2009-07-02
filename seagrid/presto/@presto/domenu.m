function theResult = domenu(self, theEvent, varargin)

% presto/domenu -- Handler for "presto" menus.
%  domenu(self, 'theEvent') handles 'theEvent' on
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

try
	theHandler = handler(self, theEvent);
	result = builtin('feval', theHandler, self, theEvent);
catch
	errormsg(self)
end

if nargout > 0
	theResult = self;
end
