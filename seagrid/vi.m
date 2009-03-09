function vi(theFile)

% vi -- Invoke "vi" on Unix.
%  vi('theFile') invokes "vi" on 'theFile' in Unix.
%   If 'theFile' is empty or a wildcard, the "uigetfile"
%   dialog is activated to get the actual filename.  If
%   the OS is not Unix, "edit" is invoked.
 
% svn $Id$
%=======================================================================
% Copyright (C) 1999 Dr. Charles R. Denham, ZYDECO.
%  All Rights Reserved.
%   Disclosure without explicit written consent from the
%    copyright owner does not constitute publication.
%=======================================================================
 
% Version of 01-Jun-1999 09:04:48.
% Updated    06-Mar-2000 10:50:08.

if nargin < 1, theFile = '*'; end

if any(theFile == '*') | isempty(theFile)
	[f, p] = uigetfile(theFile);
	if ~any(f), return, end
	if p(end) ~= filesep, p(end+1) = filesep; end
	theFile = [p f];
end

w = which(theFile);

if any(w)
	clear(w)
	if isunix
		eval(['!vi ' w])
	else
		edit(w)
	end
else
	disp([' ## No such file: ' theFile])
end
