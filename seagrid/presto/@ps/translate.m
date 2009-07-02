function theResult = translate(self, theEvent)

% ps/translate -- Translate a "ps" event.
%  translate(self, 'theEvent') converts 'theEvent'
%   to canonical form for use by "ps" event
%   handlers.  The result is lowercase, free of
%   blanks, and consists entirely of alphanumeric
%   characters and '_'.
 
% svn $Id$
%=======================================================================
% Copyright (C) 1999 Dr. Charles R. Denham, ZYDECO.
%  All Rights Reserved.
%   Disclosure without explicit written consent from the
%    copyright owner does not constitute publication.
%=======================================================================
 
% Version of 09-Dec-1999 02:45:49.
% Updated    09-Dec-1999 02:45:49.

if nargout > 0, theResult = []; end
if nargin < 2, help(mfilename), return, end

result = lower(theEvent);

theLetters = char([abs('a'):abs('z') abs('0'):abs('9') '_']);
for i = length(result):-1:1
	if ~any(result(i) == theLetters)
		result(i) = '';
	end
end

if nargout > 0
	theResult = result;
else
	assignin('caller', 'ans', result)
	disp(result)
end
