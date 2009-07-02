function theResult = help_handlers(self, varargin)

%    Event-handlers are registered with "ps" using the
% "handler" method, which requires the name of an event,
% such as 'Print', and the name of the corresponding method
% that will handle the event, such as 'doprint'.  Any number
% of event/handler pairs can be placed in the argument list,
% as in "self = handler(self, theEvent, theHandler, ...)".
%    Handlers are always called with three input arguments,
% as in "self = doprint(self, theEvent, theMessage)".
%    See "help ps/handler".
 
% svn $Id$
%=======================================================================
% Copyright (C) 1999 Dr. Charles R. Denham, ZYDECO.
%  All Rights Reserved.
%   Disclosure without explicit written consent from the
%    copyright owner does not constitute publication.
%=======================================================================
 
% Version of 04-Nov-1999 14:41:34.
% Updated    05-Nov-1999 20:08:50.

h = help(mfilename);
helpdlg(h, 'PS Handlers')

if nargout > 0, theResult = self; end
