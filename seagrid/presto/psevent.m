function theResult = psevent(theEvent, theMessage)

% psevent -- Gateway for "ps" events.
%  psevent('theEvent', theMessage) calls the "doevent"
%   method of the "ps" object associated with the gcbf
%   (current-callback-figure), if any.
 
% svn $Id$
%=======================================================================
% Copyright (C) 1999 Dr. Charles R. Denham, ZYDECO.
%  All Rights Reserved.
%   Disclosure without explicit written consent from the
%    copyright owner does not constitute publication.
%=======================================================================
 
% Version of 09-Dec-1999 02:31:41.
% Updated    09-Dec-1999 02:31:41.

if nargout > 0, theResult = []; end
if nargin < 1, theEvent = 'Callback'; end
if nargin < 2, theMessage = []; end

if isempty(gcbf), return, end

if nargout > 0
	theResult = doevent(ps(gcbf), theEvent, theMessage);
else
	doevent(ps(gcbf), theEvent, theMessage);
end
