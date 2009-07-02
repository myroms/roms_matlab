function theResult = putfile(self, theFileFilter, thePrompt)

% presto/putfile -- Invoke "uiputfile".
%  putfile(self, 'theSuggested', thePrompt) invokes the
%   "uiputfile" dialog for the given suggested filename
%   and prompt.

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
if nargin < 2, theSuggested = 'unnamed.mat'; end
if nargin < 3, thePrompt = 'Select Output File'; end

[theFile, thePath] = uiputfile(theSuggested, thePrompt);
if ~any(theFile), return, end

if thePath(end) ~= filesep, thePath(end+1) = filesep; end
result = [thePath theFile];

if nargout > 0
	theResult = result;
else
	assignin('caller', 'ans', result)
	disp(result)
end
