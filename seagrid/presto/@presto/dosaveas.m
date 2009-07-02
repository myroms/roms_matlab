function theResult = dosaveas(self, varargin)

% presto/dosaveas -- Save "presto" data.
%  dosaveas(self) saves the data of self,
%   a "presto" object, to the chosen file.

% svn $Id$
%=======================================================================
% Copyright (C) 1999 Dr. Charles R. Denham, ZYDECO.
%  All Rights Reserved.
%   Disclosure without explicit written consent from the
%    copyright owner does not constitute publication.
%=======================================================================
 
% Version of 02-Nov-1999 23:17:00.
% Updated    02-Nov-1999 23:17:00.

if nargout > 0, theResult = []; end
if nargin < 1, help(mfilename), return, end

self = dosave(self, '*');

if nargout > 0
	theResult = self;
end
