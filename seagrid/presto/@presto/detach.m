function theResult = detach(self)

% presto/detach -- Remove "presto" object from handle.
%  detach(self) removes self, a "presto" object, from
%   its associated handle.

% svn $Id$
%=======================================================================
% Copyright (C) 1999 Dr. Charles R. Denham, ZYDECO.
%  All Rights Reserved.
%   Disclosure without explicit written consent from the
%    copyright owner does not constitute publication.
%=======================================================================
 
% Version of 01-Nov-1999 10:45:28.
% Updated    01-Nov-1999 10:45:28.

if nargout > 0, theResult = presto(self); end
if nargin < 1, help(mfilename), return, end

set(handle(self), 'UserData', self.UserData)
