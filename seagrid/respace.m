function theResult = respace(x, x1, x2)

% respace -- Re-space values linearly.
%  respace(x, x1, x2) respaces the values of vector x,
%   such that the original x1 is moved to x2 and all
%   the others are correspondingly moved linearly
%   between the end-points.
 
% svn $Id$
%=======================================================================
% Copyright (C) 1999 Dr. Charles R. Denham, ZYDECO.
%  All Rights Reserved.
%   Disclosure without explicit written consent from the
%    copyright owner does not constitute publication.
%=======================================================================
 
% Version of 23-Apr-1999 13:21:29.
% Updated    25-Aug-1999 07:16:00.

if nargin < 3, help(mfilename), return, end

xmin = min(x(:)); x = x - xmin;
xmax = max(x(:)); x = x ./ xmax;

x1 = (x1 - xmin); x1 = x1 ./ xmax;
x2 = (x2 - xmin); x2 = x2 ./ xmax;

i = find(x <= x1);
j = find(x > x1);

x(i) = x(i) .* x2 ./ x1;
x(j) = 1 - (1 - x(j)) .* (1 - x2) ./ (1 - x1);

result = x .* xmax + xmin;

if nargout > 0
	theResult = result;
else
	disp(result)
end
