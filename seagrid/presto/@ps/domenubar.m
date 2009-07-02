function theResult = domenubar(self, varargin)

% ps/domenubar -- Toggle the "ps" menubar.
%  domenubar(self) toggles the menubar on behalf
%   of self, a "ps" object.
 
% svn $Id$
%=======================================================================
% Copyright (C) 1999 Dr. Charles R. Denham, ZYDECO.
%  All Rights Reserved.
%   Disclosure without explicit written consent from the
%    copyright owner does not constitute publication.
%=======================================================================
 
% Version of 04-Nov-1999 15:32:53.
% Updated    04-Nov-1999 15:32:53.

if nargout > 0, theResult = []; end
if nargin < 1, help(mfilename), return, end

theMenu = gcbo;
theFigure = gcbf;

switch lower(get(theFigure, 'MenuBar'))
case 'figure'
	set(theFigure, 'MenuBar', 'none')
	if any(theMenu), set(theMenu, 'Checked', 'off'); end
case 'none'
	set(theFigure, 'MenuBar', 'figure')
	if any(theMenu), set(theMenu, 'Checked', 'on'); end
end

if nargout > 0, theResult = self; end
