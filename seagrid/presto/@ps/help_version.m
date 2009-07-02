function theResult = help_version(self, varargin)

% presto/help_version-- Show "presto" version.
%  help_version(self) shows the version of 'presto"
%   currently active.

% svn $Id$
%=======================================================================
% Copyright (C) 1999 Dr. Charles R. Denham, ZYDECO.
%  All Rights Reserved.
%   Disclosure without explicit written consent from the
%    copyright owner does not constitute publication.
%=======================================================================
 
% Version of 04-Nov-1999 15:17:29.
% Updated    04-Nov-1999 15:17:29.

h = help('ps/version');
helpdlg(h, 'PS Version')

if nargout > 0, theResult = self; end
