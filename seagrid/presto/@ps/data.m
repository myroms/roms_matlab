function theResult = data(self, theNewData)

% ps/data -- Set/get "ps" data.
%  data(self) gets the data associated
%   with self, a "ps" object.
%  data(self, theNewData) sets the data
%   associated with self to theNewData.
 
% svn $Id$
%=======================================================================
% Copyright (C) 1999 Dr. Charles R. Denham, ZYDECO.
%  All Rights Reserved.
%   Disclosure without explicit written consent from the
%    copyright owner does not constitute publication.
%=======================================================================
 
% Version of 09-Dec-1999 00:34:18.
% Updated    09-Dec-1999 05:26:41.

if nargout > 0, theResult = []; end
if nargin < 1, help(mfilename), return, end

theHandle = handle(self);

u = get(theHandle, 'UserData');

if nargin < 2
	result = u.ps_Data;
else
	u.ps_Data = theNewData;
	set(theHandle, 'UserData', u)
	result = self;
end

if nargin > 0
	theResult = result;
else
	assignin('caller', 'ans', result)
	disp(result)
end
