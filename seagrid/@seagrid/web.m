function theResult = web(self)

% web -- World Wide Web site of the SeaGrid Toolbox.
%  web(no argument) displays or returns the WWW
%   site for the SeaGrid Toolbox.  If displayed,
%   a dialog asks whether to go there.
 
% svn $Id$
%=======================================================================
% Copyright (C) 1999 Dr. Charles R. Denham, ZYDECO.
%  All Rights Reserved.
%   Disclosure without explicit written consent from the
%    copyright owner does not constitute publication.
%=======================================================================
 
% Version of 28-Apr-1999 19:30:28.
% Updated    29-Sep-1999 11:19:35.

theURL = 'http://crusty.er.usgs.gov/~cdenham/seagrid/seagrid.html';

if nargout > 0
	theResult = theURL;
else
	disp(['## SeaGrid Toolbox Home Page:'])
	disp(['## ' theURL])
	theButton =  questdlg('Go To SeaGrid Toolbox Home Page?', 'WWW', 'Yes', 'No', 'No');
	if isequal(theButton, 'Yes')
		theStatus = web(theURL);
		switch theStatus
		case 1
			disp(' ## Could not find Web Browser.')
			disp(' ## See "help web".')
		case 2
			disp(' ## Web Browser found, but could not be launched.')
			disp(' ## See "help web".')
			help(web)
		otherwise
		end
	end
end
