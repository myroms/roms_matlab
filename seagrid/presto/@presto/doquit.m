function theResult = doquit(self, varargin)

% presto/doquit -- Quit a "proxy" process.
%  doquit(self) deletes the handle associated
%   with self, a "presto" object, then returns
%   the empty-matrix [].

% svn $Id$
%=======================================================================
% Copyright (C) 1999 Dr. Charles R. Denham, ZYDECO.
%  All Rights Reserved.
%   Disclosure without explicit written consent from the
%    copyright owner does not constitute publication.
%=======================================================================
 
% Version of 02-Nov-1999 23:18:31.
% Updated    02-Nov-1999 23:18:31.

if nargout > 0, theResult = []; end
if nargin < 1, help(mfilename), return, end

delete(gcbf)

if nargout > 0
	theResult = [];
end
