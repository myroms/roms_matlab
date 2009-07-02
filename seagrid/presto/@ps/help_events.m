function theResult = help_events(self, varargin)

%    Whenever a mouse event occurs on a "ps" object,
% a callback reading "psevent theEventName" is executed,
% where theEventName will be something like "ButtonDownFcn".
% The "ps/doevent" method is responsible for translating
% the event into an action.  Events from controls and menus
% are converted to the name of the callback-object.  All
% events are then converted to lower-case, deblanked, and
% trimmed, leaving only alphanumeric characters.  The
% result is given to a 'switch" ladder for final dispatching.
%    When a new class is derived from "ps", it must
% supply a "doevent" method of its own to override the
% default "ps" behavior, as appropriate.  To invoke
% the default behavior in a particular situation, call
% "inherit(self, 'theEventName')".
%    See "help psevent" and "ps/doevent".
 
% svn $Id$
%=======================================================================
% Copyright (C) 1999 Dr. Charles R. Denham, ZYDECO.
%  All Rights Reserved.
%   Disclosure without explicit written consent from the
%    copyright owner does not constitute publication.
%=======================================================================
 
% Version of 04-Nov-1999 14:41:34.
% Updated    04-Nov-1999 14:54:02.

h = help(mfilename);
helpdlg(h, 'PS Events')

if nargout > 0, theResult = self; end
