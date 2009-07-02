function [theResult, theEvent] = handler(self, varargin)

% ps/handler -- Register a "ps" event-handler.
%  handler(self, 'theEvent', 'theHandler') registers 'theEvent'
%   and 'theHandler' on behalf of self, a "ps" object.
%   Additional event/handler pairs can be given in the
%   argument-list.
%  handler(self, 'theEvent') returns the handler for theEvent,
%   or [] is no such handler has been registered.

% svn $Id$
%=======================================================================
% Copyright (C) 1999 Dr. Charles R. Denham, ZYDECO.
%  All Rights Reserved.
%   Disclosure without explicit written consent from the
%    copyright owner does not constitute publication.
%=======================================================================
 
% Version of 05-Nov-1999 00:21:39.
% Updated    05-Nov-1999 00:21:39.

if nargout > 0, theResult = []; end
if nargin < 1, help(mfilename), return, end

theHandlers = psget(self, 'itsHandlers');

% Get all the handlers.

if nargin < 2
	if nargout > 0
		theResult = theHandlers;
	else
		assignin('caller', 'ans', theHandlers)
		disp(theHandlers)
	end
	return
end

% Clean-up the event string.

for k = 1:2:length(varargin)
	varargin{k} = translate(self, varargin{k});
end

% Return the handler for the event.

if nargin < 3
	theEvent = varargin{1};
	theHandler = [];
	if ~isempty(theHandlers)
		if isfield(theHandlers, theEvent)
			theHandler = getfield(theHandlers, theEvent);
		end
	end
	if nargout > 0
		theResult = theHandler;
	else
		assignin('caller', 'ans', theHandler)
		disp(theHandler)
	end
	return
end

% Register the events and handlers.

for k = 1:2:length(varargin)
	theEvent = varargin{k};
	theHandler = varargin{k+1};
	theHandlers = setfield(theHandlers, theEvent, theHandler);
end

self = psset(self, 'itsHandlers',  theHandlers);

if nargout > 0
	theResult = self;
else
	assignin('caller', 'ans', self)
end
