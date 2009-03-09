function theResult = rescale(x, xmin, xmax)

% rescale -- Scale data to a range.
%  rescale(x, xmin, xmax) rescales the values of x
%   to the range [xmin:xmax] (default = [0:1]).
%   (Needs upgrade to handle constant data.)
 
% svn $Id$
%=======================================================================
% Copyright (C) 1999 Dr. Charles R. Denham, ZYDECO.
%  All Rights Reserved.
%   Disclosure without explicit written consent from the
%    copyright owner does not constitute publication.
%=======================================================================
 
% Version of 23-Apr-1999 13:19:07.

if nargin < 1, help(mfilename), return, end

if nargin < 2, xmin = 0; end
if nargin < 3, xmax = 1; end

x(:) = x - min(x(:));
x(:) = x ./ max(x(:));

x = x .* (xmax - xmin) + xmin;

if nargout > 0
	theResult = x;
else
	disp(x)
end
