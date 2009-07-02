function theResult = dorevert(self, varargin)

% presto/dorevert -- Revert to saved file.
%  dorevert(self) reverts to the saved file
%   to reset the contents of self.

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

theFile = '';
if isfield(self, 'itsSavedFile')
	theFile = self.itsSavedFile;
end

if any(theFile)
	self = doopen(self, theFile);
end

if nargout > 0
	theResult = self;
end
