function theResult = errormsg(self, theMessage)

% presto/errormsg -- Error message for "presto".
%  errormsg(self, 'theMessage') posts a message-box
%   with the given message (default = lasterr).
 
% svn $Id$
%=======================================================================
% Copyright (C) 1999 Dr. Charles R. Denham, ZYDECO.
%  All Rights Reserved.
%   Disclosure without explicit written consent from the
%    copyright owner does not constitute publication.
%=======================================================================
 
% Version of 05-Nov-1999 09:45:13.
% Updated    05-Nov-1999 09:45:13.

PAUSE_TIME = 3;

if nargout > 0, theResult = []; end
if nargin < 1, help(mfilename), return, end
if nargin < 2, theMessage = lasterr; end

theMessageBox = msgbox(theMessage, 'warn');
pause(PAUSE_TIME)
if ishandle(theMessageBox), delete(theMessageBox), end

if nargout > 0
	theResult = self;
end
