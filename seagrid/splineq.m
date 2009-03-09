function [ti, pp] = splineq(t, z, nPoints, theSlopeFlag)

% splineq -- Equally-spaced spline interpolation.
%  [ti, pp] = splineq(t, z, nPoints, theSlopeFlag) interpolates
%   complex z(t) at nPoints equally spaced along the track of
%   the spline curve, including the first and last points.  The
%   returned ti and pp are such that ppval(pp, ti) will provide
%   the actual equally-spaced values.  If theSlopeFlag is TRUE,
%   the end-slopes of the spline are pinned to the linear trend
%   between the first two and the last two points, respectively.
%  zi = splineq(t, z, nPoints, theSlopeFlag) returns the actual
%  interpolated values.
%  [ti, pp] = splineq(t, z, ti, theSlopeFlag) interpolates the
%  z(t) curve at the relative positions [0...1] along the
%  length of the curve, given by ti.
% zi = splineq(t, z, ti, theSlopeFlag) returns the actual
%  interpolated values.
% splineq(m, n) demonstrates itself with m random data points,
%  interpolated at n positions (default = 10*m), using pinned
%  slopes.  Click on the curve to do another demonstration.
%  The "zoomsafe" facility is invoked if available.
% spline('demo') calls "spline(5, 50)".
 
% svn $Id$
%=======================================================================
% Copyright (C) 1999 Dr. Charles R. Denham, ZYDECO.
%  All Rights Reserved.
%   Disclosure without explicit written consent from the
%    copyright owner does not constitute publication.
%=======================================================================
 
% Version of 26-May-1999 15:37:29.
% Updated    31-Aug-1999 13:45:25.

if nargin < 1, help(mfilename), t = 'demo'; end

if isequal(t, 'demo'), t = 5; end
if ischar(t), t = eval(t); end

if length(t) == 1
	m = t;
	n = 10*m;
	if nargin > 1, n = z; end
	if ischar(n), n = eval(n); end
	x = sort(rand(1, m));
	x = x - min(x);
	x = x / max(x);	y = rand(size(x));
	z = x + sqrt(-1)*y;
	t = linspace(min(x), max(x), m);
	flag = 1;
	if (1)
		[tti, pp] = feval(mfilename, t, z, n, flag);
	else   % Random interpolates.
		r = sort(rand(1, n));
		r(1) = 0; r(end) = 1;
		[tti, pp] = feval(mfilename, t, z, r, flag);
	end
	zzi = ppval(pp, tti);
	xxi = real(zzi); yyi = imag(zzi);
	h = plot(x, y, 'r*', xxi, yyi, 'bo-');
	legend(h, 'data', 'interp')
	theBDF = [mfilename '(' int2str(m) ', ' int2str(n) ')'];
	set(h, 'ButtonDownFcn', theBDF)
	title(theBDF), xlabel('x'), ylabel('y')
	axis equal
	eval('zoomsafe', ' ')
	figure(gcf)
	set(gcf, 'Name', 'SplineQ Demo')
	return
end

if nargin < 2, help(mfilename), return, end
if nargin < 3, nPoints = 10 * length(z); end
if nargin < 4, theSlopeFlag = 0; end

[oldm, n] = size(z);
if oldm == 1, z = z(:); end

t = t(:);

multiplier = 10;

len = length(nPoints);
if len < 2
	ttemp = linspace(min(t), max(t), multiplier*nPoints+1).';
	di = linspace(0, 1, nPoints);
else
	ttemp = linspace(min(t), max(t), multiplier*len+1).';
	di = nPoints;
end

if theSlopeFlag
	slope = diff(z) ./ diff(t);
	z = [slope(1); z; slope(end)];
end

pp = spline(t, z);
ztemp = ppval(pp, ttemp);

d = [0; cumsum(abs(diff(ztemp)))];
d = d - min(d);
d = d / max(d);

ti = spline(d, ttemp, di);
if oldm == 1, ti = ti.'; end

if nargout < 2
	ti = ppval(pp, ti);
end
