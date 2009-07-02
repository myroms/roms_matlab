function theResult = doevent(self, theEvent, theMessage, varargin)

% ps/doevent -- Call "ps" event handler.
%  doevent(self, theEvent, theMessage) calls the registered
%   event handler on behalf of self, a "ps" object.  Menus
%   and controls are processed according to their "Tag"
%   property.  A notice is posted whenever an appropriate
%   event-handler cannot be found.
 
% svn $Id$
%=======================================================================
% Copyright (C) 1999 Dr. Charles R. Denham, ZYDECO.
%  All Rights Reserved.
%   Disclosure without explicit written consent from the
%    copyright owner does not constitute publication.
%=======================================================================
 
% Version of 05-Nov-1999 19:28:46.
% Updated    17-Dec-1999 16:56:40.

if nargout > 0, theResult = []; end
if nargin < 1, help(mfilename), return, end
if nargin < 2, theEvent = '(none)'; end
if nargin < 3, theMessage = []; end

result = [];

h = gcbo;
if isempty(h) & ishandle(theMessage)
	h = theMessage;
end

% Menus and controls are processed via "Tag".

switch get(h, 'Type')
case {'uimenu', 'uicontrol'}
	theEvent = get(h, 'Tag');
end

theHandler = handler(self, theEvent);

okay = logical(0);

if ~isempty(theHandler)
	try
		result = builtin('feval', theHandler, self, theEvent, theMessage);
		okay = logical(1);
	catch
		disp([' ## Error while executing handler: ' theHandler])
	end
end

if ~okay
	theEvent = translate(self, theEvent);
	if isempty(theEvent), theEvent = '(none)'; end
	theType = get(gcbo, 'Type');
	disp([' ## No valid event-handler for "' class(self) '" class: "' ...
				theEvent '" (' theType ' = ' num2str(gcbo) ')'])
end

if nargout > 0, theResult = self; end
