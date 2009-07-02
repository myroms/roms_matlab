function theResult = event(self, theEvent, theMessage)

% pst/event -- Event dispatcher for "presto".
%  event(self, theEvent, theMessage) dispatches theEvent
%   to the appropriate handler on behalf of self, a "pst"
%   object.  In its present form, this routine only calls
%   event-handlers inherited from "presto", its superclass.
%   Users must modify this method to customize its behavior.
 
% svn $Id$
%=======================================================================
% Copyright (C) 1999 Dr. Charles R. Denham, ZYDECO.
%  All Rights Reserved.
%   Disclosure without explicit written consent from the
%    copyright owner does not constitute publication.
%=======================================================================
 
% Version of 05-Nov-1999 20:11:19.
% Updated    05-Nov-1999 21:04:58.

if nargout > 0, theResult = []; end
if nargin < 1, help(mfilename), return, end
if nargin < 2, theEvent = ''; end
if nargin < 3, theMessage = ''; end

result = [];

% Get a genuine copy of self.

self = retrieve(self);

% Dispatch the event.

theEvent = translate(self, theEvent);

switch theEvent
case {'some custom event'}   % Customized processing here.
	result = doevent(self, theEvent, theMessage);
otherwise   % Inherit.
	result = inherit('event', self, theEvent);
end

% Restore self.

if isa(result, 'pst')
	self = result;
	self = restore(self);
end

% Output.

if nargout > 0, theResult = self; end
