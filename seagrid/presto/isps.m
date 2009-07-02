function theResult = isps(theItem)

% isps -- Is this a "ps" object?
%  isps(theItem) returns logical(1) (TRUE) if theItem
%   is a "ps" object or is derived from the "ps"
%   class; otherwise, it returns logical(0) (FALSE).
%   If theItem is a handle, the corresponding "UserData"
%   is examined instead.

% svn $Id$
%=======================================================================
% Copyright (C) 1999 Dr. Charles R. Denham, ZYDECO.
%  All Rights Reserved.
%   Disclosure without explicit written consent from the
%    copyright owner does not constitute publication.
%=======================================================================
 
% Version of 07-Dec-1999 23:19:23.
% Updated    09-Dec-1999 05:26:41.

if nargout > 0, theResult = []; end
if nargin < 1, help(mfilename), return, end

result = logical(0);

switch class(theItem)
case 'double'
	if ishandle(theItem)
		u = get(theItem, 'UserData');
		if isstruct(u) & isfield(u, 'ps_Self')
			result = isa(getfield(u, 'ps_Self'), 'ps');
		end
	end
otherwise
	result = isa(theItem, 'ps');
end

if nargout > 0
	theResult = result;
else
	assignin('caller', 'ans', result)
	disp(result)
end
