function theResult = psbind(self)

% psbind -- Bind a "ps" object to its handle.
%  psbind(self) binds self, a "ps" object, to
%   its handle.
 
% svn $Id$
%=======================================================================
% Copyright (C) 1999 Dr. Charles R. Denham, ZYDECO.
%  All Rights Reserved.
%   Disclosure without explicit written consent from the
%    copyright owner does not constitute publication.
%=======================================================================
 
% Version of 10-Dec-1999 15:03:28.
% Updated    10-Dec-1999 15:03:28.

if nargout > 0, theResult = []; end
if nargin < 1, help(mfilename), return, end

theHandle = handle(self);

u = [];
u.ps_Self = self;
u.ps_Data = [];
u.ps_Data.Handle = theHandle;
u.ps_Data.UserData = get(theHandle, 'UserData');

set(theHandle, 'UserData', u)

if nargout > 0
	theResult = self;
else
	assignin('caller', 'ans', self)
end
