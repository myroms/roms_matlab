function theResult = doscroll(self, theEvent, theMessage)

% ps/doscroll -- Handler for "ps" scrollbars.
%  doscroll(self, theEvent, theMessage) pans the axes
%   in keeping with the active scrollbars, on behalf
%   of self, a "ps" object.

% svn $Id$
%=======================================================================
% Copyright (C) 1999 Dr. Charles R. Denham, ZYDECO.
%  All Rights Reserved.
%   Disclosure without explicit written consent from the
%    copyright owner does not constitute publication.
%=======================================================================
 
% Version of 05-Nov-1999 21:34:53.
% Updated    05-Nov-1999 21:34:53.

if nargout > 0, theResult = []; end
if nargin < 1, help(mfilename), return, end

theXLim = get(gca, 'XLim');   % bottom scrollbar.
theYLim = get(gca, 'YLim');   % right scrollbar.
theZLim = get(gca, 'ZLim');   % left scrollbar.
theCLim = get(gca, 'CLim');   % top scrollbar.

axis tight
set(gca, 'CLimMode', 'auto')

theEvent = translate(self, theEvent);

switch theEvent
case 'bottom'
	theXLim = get(gca, 'XLim');
	smin = get(gcbo, 'Min');
	smax = get(gcbo, 'Max');
	value = get(gco, 'Value');
	frac = value / (smax - smin);
	theXLim = theXLim + (frac - 0.5)*diff(theXLim);
case 'right'
	theYLim = get(gca, 'YLim');
	smin = get(gcbo, 'Min');
	smax = get(gcbo, 'Max');
	value = get(gco, 'Value');
	frac = value / (smax - smin);
	theYLim = theYLim + (frac - 0.5)*diff(theYLim);
case 'left'
	theZLim = get(gca, 'ZLim');
	smin = get(gcbo, 'Min');
	smax = get(gcbo, 'Max');
	value = get(gco, 'Value');
	frac = value / (smax - smin);
	theZLim = theZLim + (frac - 0.5)*diff(theZLim);
case 'top'
	theCLim = get(gca, 'CLim');
	smin = get(gcbo, 'Min');
	smax = get(gcbo, 'Max');
	value = get(gco, 'Value');
	frac = value / (smax - smin);
	theCLim = theCLim + (frac - 0.5)*diff(theCLim);
end

set(gca, 'XLim', theXLim, 'YLim', theYLim, 'ZLim', theZLim, 'CLim', theCLim)

if nargout > 0
	theResult = self;
else
	assignin('caller', 'ans', self)
end
