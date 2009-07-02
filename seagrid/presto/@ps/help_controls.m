function theResult = help_controls(self, varargin)

%    The "ps/control" method adds controls to the
% "ps" window.  Standard scrollbars are indicated
% by the names: 'bottom', 'right', 'top,' or 'left'.
% Other controls are denoted by their "Style" property,
% as in 'pushbutton'.
%    Each control requires a "layout", consisting of a
% normalized-position and a pixel-offset.  The "layout"
% information is used by "ps/doresize" to adjust
% each control whenever the window is resized.
%    The "ps/control" method returns the handle of
% the new control, so that further embellishments can
% be made.
%    See "help ps/control".
 
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
helpdlg(h, 'PS Controls')

if nargout > 0, theResult = self; end
