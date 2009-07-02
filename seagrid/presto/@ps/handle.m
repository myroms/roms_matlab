function theResult = handle(self)

% ps/handle -- Handle of a "ps" object.
%  handle(self) returns the handle associated
%   with self, a "ps" object.

% svn $Id$
%=======================================================================
% Copyright (C) 1999 Dr. Charles R. Denham, ZYDECO.
%  All Rights Reserved.
%   Disclosure without explicit written consent from the
%    copyright owner does not constitute publication.
%=======================================================================
 
% Version of 09-Dec-1999 00:30:47.
% Updated    09-Dec-1999 00:30:47.

if nargout > 0, theResult = []; end
if nargin < 1, help(mfilename), return, end

result = ps(self);

if nargout > 0
	theResult = result;
else
	assignin('caller', 'ans', result')
	disp(result)
end
