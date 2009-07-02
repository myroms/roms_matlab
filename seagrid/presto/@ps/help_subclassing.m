function theResult = help_subclassing(self, varargin)

%    Users of PS are expected to derive their own class
% from "ps".   One such class, "pst" has already been
% built as an example.  It can be copied and renamed to serve
% as a framework for the derived class.  Unless "pst/doevent"
% is modified, all events will passed intact to "ps".
%    See "help pst".

% svn $Id$
%=======================================================================
% Copyright (C) 1999 Dr. Charles R. Denham, ZYDECO.
%  All Rights Reserved.
%   Disclosure without explicit written consent from the
%    copyright owner does not constitute publication.
%=======================================================================
 
% Version of 05-Nov-1999 20:46:35.
% Updated    05-Nov-1999 20:46:35.

h = help(mfilename);
helpdlg(h, 'PS Subclassing')

if nargout > 0, theResult = self; end
