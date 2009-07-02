function theResult = doquit(self, varargin)

% ps/doquit -- Quit a "ps" process.
%  doquit(self) deletes the handle associated
%   with self, a "ps" object, then returns
%   the empty-matrix [].
 
% svn $Id$
%=======================================================================
% Copyright (C) 1999 Dr. Charles R. Denham, ZYDECO.
%  All Rights Reserved.
%   Disclosure without explicit written consent from the
%    copyright owner does not constitute publication.
%=======================================================================
 
% Version of 02-Nov-1999 23:18:31.
% Updated    09-Dec-1999 03:13:42.

if nargout > 0, theResult = []; end
if nargin < 1, help(mfilename), return, end

delete(gcbf)

if nargout > 0
	theResult = [];
end
