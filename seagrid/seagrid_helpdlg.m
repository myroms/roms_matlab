function handle = seagrid_helpdlg(helpstring, dlgname, color)

% seagrid_helpdlg -- "helpdlg" + background color.
%  seagrid_helpdlg(helpstring, dlgname, color) calls "helpdlg",
%   then adjusts the background to the given color.
 
% svn $Id$
%=======================================================================
% Copyright (C) 1999 Dr. Charles R. Denham, ZYDECO.
%  All Rights Reserved.
%   Disclosure without explicit written consent from the
%    copyright owner does not constitute publication.
%=======================================================================
 
% Version of 24-Aug-1999 14:57:37.
% Updated    24-Aug-1999 14:57:37.

if nargin < 1
	help(mfilename)
	h = seagrid_helpdlg(help(mfilename));
	if nargout > 0, handle = h; end
	return
end
if nargin < 2, dlgname = 'SeaGrid Help'; end
if nargin < 3, color = [1 1 1]; end

h = helpdlg(helpstring, dlgname);
set(h, 'Color', color)

t = findobj(h, 'Type', 'uicontrol');
set(t, 'BackgroundColor', color)

if nargout > 0, handle = h; end
