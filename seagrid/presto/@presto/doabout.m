function doabout(self)

%    The "presto" class provides a framework for
% processing interactive graphical events.  It
% features "dot" syntax for the setting/getting
% of graphical and user-defined properties.  All
% events are trapped by calls to "event", whose
% sole argument is the name of Matlab callback
% (e.g. "ButtonDownFcn").  Use "menu" to add
% menu items, "control" to add controls, and
% "enable" to automatically enable the relevant
% callbacks.  The file "presto_test.m" gives a
% demonstration.
%    The "presto" class should be super-classed,
% and the "event" method should be overloaded.

% svn $Id$
%=======================================================================
% Copyright (C) 1999 Dr. Charles R. Denham, ZYDECO.
%  All Rights Reserved.
%   Disclosure without explicit written consent from the
%    copyright owner does not constitute publication.
%=======================================================================
 
% Version of 01-Nov-1999 13:53:56.
% Updated    04-Nov-1999 15:56:44.

helpdlg(help(mfilename), 'About Presto')
