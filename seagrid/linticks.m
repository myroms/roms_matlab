function theResult = linticks(theMin, theMax, theTickCount)

% linticks -- Tick positions within a range.
%  linticks(theMin, theMax, theTickCount) returns a vector
%   of tick positions that span the interval from theMin
%   to theMax, with theTickCount or fewer elements.
%   Based on the "goodscales" scheme used by the Matlab
%   "plotyy" routine.
 
% svn $Id$
%=======================================================================
% Copyright (C) 2000 Dr. Charles R. Denham, ZYDECO.
%  All Rights Reserved.
%   Disclosure without explicit written consent from the
%    copyright owner does not constitute publication.
%=======================================================================
 
% Version of 08-Feb-2000 21:27:14.
% Updated    08-Feb-2000 21:27:14.

if nargout > 0, theResult = []; end

if nargin < 1, help(mfilename), return, end

if ischar(theMin), theMin = eval(theMin); end
if ischar(theMax), theMax = eval(theMax); end
if ischar(theTickCount), theTickCount = eval(theTickCount); end

[low, high, ticks] = goodscales(theMin, theMax);

f = find(ticks <= theTickCount);
if any(f)
	f = f(1);
	result = linspace(low(f), high(f), ticks(f));
	result(abs(result) < sqrt(eps)) = 0;
end

if nargout > 0
	theResult = result;
else
	assignin('caller', 'ans', result)
	disp(result)
end


% ---------- goodscales ---------- %


function [low, high, ticks] = goodscales(xmin, xmax)

% GOODSCALES -- Returns parameters for "good" scales.
%  [LOW, HIGH, TICKS] = GOODSCALES(XMIN, XMAX) returns lower and upper
%   axis limits (LOW and HIGH) that span the interval (XMIN, XMAX) 
%   with "nice" tick spacing.  The number of major axis ticks is 
%   returned in TICKS.

% ** Liberated from "plotyy.m".

% Modifications:
% Copyright (C) 2000 Dr. Charles R. Denham, ZYDECO.
%  All Rights Reserved.
%   Disclosure without explicit written consent from the
%    copyright owner does not constitute publication.
 
% Version of 08-Feb-2000 20:50:36.
% Updated    08-Feb-2000 21:16:45.

if nargin < 1, help(mfilename), return, end

if ischar(xmin), xmin = eval(xmin); end
if ischar(xmax), xmax = eval(xmax); end

BestDelta = [ 0.1 0.2 0.5 1 2 5 10 20 50 ];

xmin = min(xmin(:));
xmax = max(xmax(:));

if xmin == xmax
	lo = xmin;
	high = xmax + 1;
	ticks = 1;
else
	Xdelta = xmax-xmin;
	delta = 10.^(round(log10(Xdelta) - 1)) * BestDelta;
	high = delta .* ceil(xmax ./ delta);
	lo = delta .* floor(xmin ./ delta);
	ticks = round((high - lo) ./ delta) + 1;
end

if nargout > 0
	low = lo;
elseif nargout == 1
	low = [lo; high; ticks];
else
	disp([lo; high; ticks])
end
