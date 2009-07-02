function theResult = doevent(self, theEvent, theMessage, varargin)

% presto/doevent -- Call "presto" event handler.
%  doevent(self, theEvent, theMessage) calls the appropriate
%   event handler on behalf of self, a "presto" object.
 
% svn $Id$
%=======================================================================
% Copyright (C) 1999 Dr. Charles R. Denham, ZYDECO.
%  All Rights Reserved.
%   Disclosure without explicit written consent from the
%    copyright owner does not constitute publication.
%=======================================================================
 
% Version of 05-Nov-1999 19:28:46.
% Updated    08-Dec-1999 11:47:11.

if nargout > 0, theResult = []; end
if nargin < 1, help(mfilename), return, end
if nargin < 2, theEvent = '(none)'; end
if nargin < 3, theMessage = []; end

result = [];

switch get(gcbo, 'Type')
case 'uimenu*'
	result = domenu(self, theEvent, theMessage);
case 'uicontrol*'
	result = docontrol(self, theEvent, theMessage);
otherwise
	theHandler = handler(self, theEvent);
	if ~isempty(theHandler)
		result = builtin('feval', theHandler, self, theEvent, theMessage);
	else
		if isempty(theEvent), theEvent = '(none)'; end
		theType = get(gcbo, 'Type');
		disp([' ## ' datestr(now, 13) ' ' theEvent ' ' theType ' ' num2str(gcbo)])
	end
end

if nargout > 0
	theResult = self;
end
