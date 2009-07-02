function theResult = dohome(self, varargin)

% ps/dohome -- Set controls to "home" values.
%  dohome(self) sets the controls associated with
%   self, a "ps" object, to their "home" value.
%   For sliders, this "Value" is midway between
%   the "Min" and "Max" properties.
%
% Note: this routine also needs to force the figure
%  to update itself.
 
% svn $Id$
%=======================================================================
% Copyright (C) 1999 Dr. Charles R. Denham, ZYDECO.
%  All Rights Reserved.
%   Disclosure without explicit written consent from the
%    copyright owner does not constitute publication.
%=======================================================================
 
% Version of 10-Dec-1999 01:08:11.
% Updated    10-Dec-1999 14:27:28.

if nargout > 0, theResult = []; end
if nargin < 1, help(mfilename), return, end

theSliders = findobj(handle(self), 'Type', 'uicontrol', ...
							'Style', 'slider');
					
for i = 1:length(theSliders)
	theMin = get(theSliders(i), 'Min');
	theMax = get(theSliders(i), 'Max');
	set(theSliders(i), 'Value', mean([theMin theMax]))
end

axis tight
set(gca, 'CLimMode', 'auto')

if nargout > 0
	theResult = self;
else
	assignin('caller', 'ans', self)
end
