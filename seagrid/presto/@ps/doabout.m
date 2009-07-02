function theResult = doabout(self, varargin)

%    The "ps" class provides a framework for
% processing interactive graphical events.  It
% features "dot" syntax for the setting/getting
% of graphical and user-defined properties.  All
% events are trapped by calls to "psevent", whose
% sole argument is the name of Matlab callback
% (e.g. "ButtonDownFcn").  Use "menu" to add
% menu items, "control" to add controls, and
% "enable" to automatically enable the relevant
% callbacks.  The file "ps_test.m" gives a
% demonstration.
%    The "ps" class should be super-classed,
% and the "event" method should be overloaded.

% ps/doabout -- Post "About PS" message.
%  doabout(self) posts the "About PS" message,
%  as seen in the block of text shown above.

% svn $Id$
%=======================================================================
% Copyright (C) 1999 Dr. Charles R. Denham, ZYDECO.
%  All Rights Reserved.
%   Disclosure without explicit written consent from the
%    copyright owner does not constitute publication.
%=======================================================================
 
% Version of 01-Nov-1999 13:53:56.
% Updated    04-Nov-1999 15:56:44.

helpdlg(help(mfilename), 'About PS')

if nargout > 0, theResult = self; end
