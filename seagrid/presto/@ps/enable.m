function theResult = enable(self)

% ps/enable -- Enable "ps" callbacks.
%  enable(self) enables all the callbacks associated
%   with self, a "ps" object, except that the "CreateFcn"
%   and "DeleteFcn" are enabled only for figures.  The
%   actions are directed to "psevent", using the actual
%   callback names, as in "psevent ButtonDownFcn".

% svn $Id$
%=======================================================================
% Copyright (C) 1999 Dr. Charles R. Denham, ZYDECO.
%  All Rights Reserved.
%   Disclosure without explicit written consent from the
%    copyright owner does not constitute publication.
%=======================================================================
 
% Version of 27-Oct-1999 23:30:10.
% Updated    09-Dec-1999 02:42:24.

if nargout > 0, theResult = []; end
if nargin < 1, help(mfilename), return, end

theHandle = handle(self);
h = findobj(theHandle);

for k = 1:length(h)
	theType = get(h(k), 'Type');
	theEvents = {};
	switch theType
	case 'figure'
		theEvents = [theEvents ...
			{'WindowButtonDownFcn', 'ResizeFcn', ...
			'CreateFcn', 'CloseRequestFcn'} ...
			];
	case {'axes', 'line', 'patch', 'surface', 'text', 'light'}
		theEvents = [theEvents {'ButtonDownFcn'}];
	case 'uicontrol'
		theEvents = [theEvents {'ButtonDownFcn', 'Callback'}];
	case 'uimenu'
		if ~any(get(h(k), 'Children'))
			theEvents = [theEvents {'Callback'}];
		end
	otherwise
	end
	for i = 1:length(theEvents)
		set(h(k), theEvents{i}, ['psevent ' theEvents{i}])
	end
end

if nargout > 0, theResult = self; end
