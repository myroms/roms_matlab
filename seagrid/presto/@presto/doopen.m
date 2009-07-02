function theResult = doopen(self, varargin)

% presto/doopen -- Open a file.
%  doopen(self) opens and reads the contents of self,
%   a "presto" object, from theFile.  If no filename
%   is given, the "uigetfile" dialog is invoked.

% svn $Id$
%=======================================================================
% Copyright (C) 1999 Dr. Charles R. Denham, ZYDECO.
%  All Rights Reserved.
%   Disclosure without explicit written consent from the
%    copyright owner does not constitute publication.
%=======================================================================
 
% Version of 01-Nov-1999 12:54:50.
% Updated    01-Nov-1999 12:54:50.

if nargout > 0, theResult = []; end
if nargin < 1, help(mfilename), return, end

theFile = '*';
if length(varargin) > 1 & ~isequal(varargin{1}, 'open')
	theFile = varargin{1};
else
	theFile = getfile(self);
end

try
	load(theFile)
	self = setps(self, 'itsSavedFile', theFile);
catch
	errormsg(self)
end

if nargout > 0
	theResult = self;
end
