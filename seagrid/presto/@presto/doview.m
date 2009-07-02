function theResult = doview(self, varargin)

% presto/view -- Toggle between view(2) and view(3).
%  view(self) toggles the view of axes associated
%   with self, a "presto" object.

% svn $Id$
%=======================================================================
% Copyright (C) 1999 Dr. Charles R. Denham, ZYDECO.
%  All Rights Reserved.
%   Disclosure without explicit written consent from the
%    copyright owner does not constitute publication.
%=======================================================================
 
% Version of 10-Nov-1999 08:17:04.
% Updated    17-Nov-1999 09:05:37.

if nargout > 0, theResult = []; end

theMenu = gcbo;

v = view;
if norm(diag(view)) > 1.9
	view(3)
	if any(theMenu), set(gcbo, 'Checked', 'on'), end
else
	view(2)
	if any(theMenu), set(gcbo, 'Checked', 'off'), end
end

if nargout > 0
	theResult = self;
else
	assignin('caller', 'ans', self)
end
