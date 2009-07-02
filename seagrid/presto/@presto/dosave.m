function theResult = dosave(self, varargin)

% presto/dosave -- Save to a file.
%  dosave(self, 'theFile') saves the contents of
%   self to theFile.  If no filename is given,
%   the "uiputfile" dialog is invoked.

% svn $Id$
%=======================================================================
% Copyright (C) 1999 Dr. Charles R. Denham, ZYDECO.
%  All Rights Reserved.
%   Disclosure without explicit written consent from the
%    copyright owner does not constitute publication.
%=======================================================================
 
% Version of 01-Nov-1999 12:54:50.
% Updated    05-Nov-1999 09:22:44.

if nargout > 0, theResult = []; end
if nargin < 1, help(mfilename), return, end

theFile = '*';
if length(varargin) > 0 & ~isequal(varargin{1}, 'save')
	theFile = varargin{1};
else
	theData = data(self);
	if isfield(theData, 'itsSavedFile')
		theFile = self.itsSavedFile;
	end
end

if any(theFile == '*')
	theFile = putfile(self);
end

try
	theData = data(self);
	theMessageBox = msgbox(['Saving to "' theFile '"']);
	save(theFile, 'theData')
	pause(1)
	if ishandle(theMessageBox), delete(theMessageBox), end
	self = setps(self, itsSavedFile, theFile);
catch
	errormsg(self)
end

if nargout > 0
	theResult = self;
end
