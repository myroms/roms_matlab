function theResult = splinesafe(t, z, ti, theSlopeFlag)

% splinesafe -- Spline interpolation with end-slope control.
%  splinesafe(t, z, theSlopeFlag) returns "pp" coefficients for
%   interpolating complex z(t), such that ppval(pp, ti) will
%   provide interpolated values at points ti.  If theSlopeFlag
%   is TRUE, the end-slopes are pinned to the linear slopes
%   between the first two and the last two points in the data,
%   respectively.  If z is a matrix, this routine works down
%   the columns, unlike "spline" itself.
%  splinesafe(t, z, ti, slopeFlag) returns the interpolates for
%   z(t) at positions ti.  The result has the same shape as ti.
%  splinesafe(nPoints) demonstrates itself with nPoints
%   (default = 5) of random data.  Click on any curve
%   to do another demonstration.  The "zoomsafe" facility
%   is active if available.  See "help zoomsafe".
%  splinesafe('demo') calls splinesafe(5).
 
% svn $Id$
%=======================================================================
% Copyright (C) 1999 Dr. Charles R. Denham, ZYDECO.
%  All Rights Reserved.
%   Disclosure without explicit written consent from the
%    copyright owner does not constitute publication.
%=======================================================================
 
% Version of 25-May-1999 08:11:28.
% Updated    31-Aug-1999 13:44:59.

if nargin < 1, help(mfilename), t = 'demo'; end
if isequal(t, 'demo'), t = 5; end
if ischar(t), t = eval(t); end

if length(t) == 1
	nPoints = t;
	t = sort(rand(1, nPoints));
	t = t - min(t); t = t / max(t);
	z = rand(size(t));
	ti = linspace(min(t), max(t), 100*length(t)+1);
	ppf = feval(mfilename, t, z, 0);   % Free-slope.
	ppp = feval(mfilename, t, z, 1);   % Pinned-slope.
	zf = ppval(ppf, ti);
	zp = ppval(ppp, ti);
	h = plot(t, z, 'bo:', ti, zf, 'r-', ti, zp, 'g-');
	theBDF = [mfilename '(' int2str(nPoints) ')'];
	set(h, 'ButtonDownFcn', theBDF)
	legend(h, 'data', 'free spline', 'pinned slope')
	title(theBDF), xlabel('t'), ylabel('z')
	figure(gcf)
	set(gcf, 'Name', 'SplineSafe Demo')
	eval('zoomsafe', ' ')
	return
end

if nargin < 3
	ti = [];
	theSlopeFlag = 0;
elseif nargin < 4
	if length(ti) == 1
		theSlopeFlag = ti;
		ti = [];
	else
		theSlopeFlag = 0;
	end
end

t = t(:);
ti = ti(:);

if min(size(z)) == 1
	z = z(:);
end

if theSlopeFlag
	dt = diff(t); dt = dt([1 end], :);
	dz = diff(z); dz = dz([1 end], :);
	[m, n] = size(dz);
	dt = dt * ones(1, n);
	slope = dz ./ dt;
	z = [slope(1, :); z; slope(end, :)];
end

% N.B. The "spline" routine works
%  across rows when z is a matrix.

pp = spline(t.', z.');   % Spline works across rows!

if isempty(ti)
	result = pp;
else
	zi = ppval(pp, ti.').';
	result = zi;
end

if nargout > 0
	theResult = result;
else
	assignin('caller', 'ans', result)
end
