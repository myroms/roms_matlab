function pp = splinep(t, z, theSlopeFlag)

% splinep -- Spline interpolation with end-slope control.
%  splinep(t, z, theSlopeFlag) returns "pp" coefficients for
%   interpolating complex z(t), such that ppval(pp, ti) will
%   provide interpolated values at points ti.  If theSlopeFlag
%   is TRUE, the end-slopes are pinned to the linear slopes
%   between the first two and the last two points in the data,
%   respectively.
%  splinep(nPoints) demonstrates itself with nPoints
%   (default = 5) of random data.  Click on any curve
%   to do another demonstration.  The "zoomsafe" facility
%   is active if available.  See "help zoomsafe".
%  splinep('demo') calls splinep(5).
 
% svn $Id$
%=======================================================================
% Copyright (C) 1999 Dr. Charles R. Denham, ZYDECO.
%  All Rights Reserved.
%   Disclosure without explicit written consent from the
%    copyright owner does not constitute publication.
%=======================================================================
 
% Version of 25-May-1999 08:11:28.
% Updated    27-May-1999 11:18:10.

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
	set(gcf, 'Name', 'SplineP Demo')
	eval('zoomsafe', ' ')
	return
end

if nargin < 3, theSlopeFlag = 0; end

if theSlopeFlag
	dt = diff(t);
	dz = diff(z);
	slope = dz ./ dt;
	z = [slope(1); z(:); slope(end)];
end

pp = spline(t(:), z(:));
