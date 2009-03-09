function showhide(theHandles)

% showhide -- Toggle visibility of a graphics object.
%  showhide(theHandles) toggles the visibility of each
%   of theHandles.
 
% svn $Id$
%=======================================================================
% Copyright (C) 1999 Dr. Charles R. Denham, ZYDECO.
%  All Rights Reserved.
%   Disclosure without explicit written consent from the
%    copyright owner does not constitute publication.
%=======================================================================
 
% Version of 21-Apr-1999 14:05:01.

if nargin < 1, help(mfilename), return, end

for i = 1:length(theHandles)
	switch get(theHandles(i), 'Visible')
	case 'on'
		theVisible = 'off';
	otherwise
		theVisible = 'on';
	end
	set(theHandles(i), 'Visible', theVisible)
end
