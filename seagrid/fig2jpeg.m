function theOutputInfo = fig2jpeg(theFilename, theFigure, theQuality)

% fig2jpeg -- Save a visible figure in JPEG format.
%  fig2jpeg('theFilename', theFigure, theQuality) saves theFigure
%   (must exist and be visible; default = gcf) to 'theFilename'
%   (default = 'unnamed.jpg'), with theQuality (default = 100).
%   If no filename is given, the Matlab "uiputfile" dialog is
%   invoked.  If an output argument is provided, the "imfinfo"
%   string is returned.
 
% svn $Id$
%=======================================================================
% Copyright (C) 2000 Dr. Charles R. Denham, ZYDECO.
%  All Rights Reserved.
%   Disclosure without explicit written consent from the
%    copyright owner does not constitute publication.
%=======================================================================
 
% Version of 03-Jan-2000 17:35:28.
% Updated    04-Jan-2000 14:03:25.

if nargout > 0, theOutputInfo = []; end

if ~any(findobj('Type', 'figure', 'Visible', 'on'))
	disp(' ## No visible figures.')
	return
end

if nargin < 1
	[f, p] = uiputfile('unnamed.jpg', 'Save Figure As JPEG:');
	if ~any(f)
		help(mfilename)
		return
	end
	if p(end) ~= filesep, p(end+1) = filesep; end
	theFilename = [p f];
end
if nargin < 2, theFigure = gcf; end
if nargin < 3, theQuality = 100; end

theOldFigure = gcf;

figure(theFigure)
[x, map] = getframe(theFigure);
imwrite(x, map, theFilename, 'jpg', 'Quality', theQuality)

figure(theOldFigure)

if nargout > 0
	theOutputInfo = imfinfo(theFilename);
end
