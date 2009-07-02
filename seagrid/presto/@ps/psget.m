function theResult = psget(self, theField)

% ps/psget -- Get a "ps" field value.
%  psget(self, 'theField') returns the value of
%   'theField' associated with self, a "ps" object.
%   The empty matrix [] is returned if no such
%   field exists.  If the field is the name of
%   a property of the handle associated with
%   self, that property is returned.
%  psget(self) returns all the current "ps" fields
%   in a "struct".
 
% svn $Id$
%=======================================================================
% Copyright (C) 1999 Dr. Charles R. Denham, ZYDECO.
%  All Rights Reserved.
%   Disclosure without explicit written consent from the
%    copyright owner does not constitute publication.
%=======================================================================
 
% Version of 09-Dec-1999 00:44:58.
% Updated    09-Dec-1999 00:44:58.

if nargout > 0, theResult = []; end
if nargin < 1, help(mfilename), return, end

result = [];

if nargin < 2
	result = psset(self);
else
	result = psset(self, theField);
end

if nargout > 0
	theResult = result;
else
	assignin('caller', 'ans', result)
end
