function theResult = getfile(self, theFileFilter, thePrompt)

% presto/getfile -- Invoke "uigetfile".
%  getfile(self, 'theFileFilter', 'thePrompt') invokes the
%   "uigetfile" dialog for the given file-filter and prompt.

% svn $Id$
%=======================================================================
% Copyright (C) 1999 Dr. Charles R. Denham, ZYDECO.
%  All Rights Reserved.
%   Disclosure without explicit written consent from the
%    copyright owner does not constitute publication.
%=======================================================================
 
% Version of 01-Nov-1999 11:37:30.
% Updated    01-Nov-1999 11:37:30.

if nargout > 0, theResult = []; end
if nargin < 1, help(mfilename), return, end
if nargin < 2, theFileFilter = '*'; end
if nargin < 3, thePrompt = 'Select Input File'; end

[theFile, thePath] = uigetfile(theFileFilter, thePrompt);
if ~any(theFile), return, end

if thePath(end) ~= filesep, thePath(end+1) = filesep; end
result = [thePath theFile];

if nargout > 0
	theResult = result;
else
	assignin('caller', 'ans', result)
	disp(result)
end
