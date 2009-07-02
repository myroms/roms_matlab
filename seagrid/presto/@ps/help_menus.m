function theResult = help_menus(self, varargin)

%    The "ps/menu" method adds menus to the PS
% menubar, using a cell-array of the desired menu labels,
% arranged in top-down fashion.  Prepended to each label
% is one ">" for each stage below the top-level menu.
% For a separator-line between menus, use "-" in place
% of one of the ">" symbols.
%    Menus call "event Callback" when selected.  Inside
% "ps/doevent", the event is translated into the menu
% label of the instigator, then processed in the "switch"
% ladder, as described in the "Help/Events" menu.
%    See "help ps/menu".

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
helpdlg(h, 'PS Menus')

if nargout < 1, theResult = self; end
