function theResult = hello(varargin)

% hello -- Display a message.
%  hello(theMessage) displays 'hello' plus theMessage.
%   If called from a subroutine, the line-number and
%   function-name are included in the display.

% svn $Id$
%=======================================================================
% Copyright (C) 1996 Dr. Charles R. Denham, ZYDECO.
%  All Rights Reserved.
%   Disclosure without written consent from the
%    copyright owner does not constitute publication.
%=======================================================================

% Updated    03-May-2001 19:11:40.

if nargin < 1, theMessage = ''; end

[stack, index] = dbstack;

if length(stack) > 1
	stack = stack(2);
else
	stack = [];
	stack.line = 0;
	stack.name = '';
end

if ~isempty(stack.name)
	s = [' ## hello at "' stack.name '" line ' int2str(stack.line)];
else
	s = [' ## hello'];
end

for i = 1:length(varargin)
   s = [s ' ' mat2str(varargin{i})];
end

if nargout < 1
	disp(s)
else
	theResult = s;
end
