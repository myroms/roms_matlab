function tick = pdf2tick(x, pdf, cdf)

% pdf2tick -- Convert probability density function to ticks.
%  pdf2tick(x, pdf, tics) returns ticks along [0..1] that
%   conform to the given pdf(x).  The number of ticks or the
%   vector of positions in the cumulative probability function
%   is set by tics (default = 101).  The pdf values represent
%   relative weights, not absolute densities.  They are fit
%   with a spline, integrated, then spline interpolated in
%   reverse to get the result.  If values are not provided
%   for x = 0 or x = 1, they are extrapolated by splines.
%  pdf2tick(ntics, npdf) demonstrates itself with ntics,
%   equally-spaced,  using npdf random pdf values.
%   Defaults are 51 and 5, respectively.
 
% svn $Id$
%=======================================================================
% Copyright (C) 1999 Dr. Charles R. Denham, ZYDECO.
%  All Rights Reserved.
%   Disclosure without explicit written consent from the
%    copyright owner does not constitute publication.
%=======================================================================
 
% Version of 13-Sep-1999 10:01:10.
% Updated    01-Oct-1999 15:36:14.

if nargin < 1, help(mfilename), x = 'demo'; end
if isequal(x, 'demo'), x = 51; end
if ischar(x), x = eval(x); end
if nargin > 1 & ischar(pdf), pdf = eval(pdf); end

isdemo = 0;

if length(x) == 1
	isdemo = 1;
	cdf = x;
	n = 5;
	if nargin > 1 & length(pdf) == 1
		n = pdf;
	end
	c = 2*n;
	x = linspace(1/c, 1-1/c, c/2);
	pdf = rand(size(x));
	pdf = pdf / max(pdf);
end

if ~isdemo & nargin < 3, cdf = 101; end

if length(cdf) == 1
	cdf = linspace(0, 1, cdf);
	cdf(1) = 0; cdf(end) = 1;
end

cdf = cdf(:).';

xi = linspace(0, 1, 5*length(cdf));   % Note 5-fold sampling.
yi = splinesafe(x, pdf, xi, 1);
scale = max(yi);

pdf = pdf / scale;   % For plotting.
yi = yi / scale;

yi(yi <= 0) = sqrt(eps);

zi = cumsum(yi);
for i = 2:length(zi)
	if zi(i) < zi(i-1)
		zi(i) = zi(i-1)
	end
end
zi = zi - min(zi);
zi = zi / max(zi);
zi(1) = 0;
zi(end) = 1;

yi = yi / scale;
pdf = pdf / scale;

t = splinesafe(zi, xi, cdf, 1);
f = find(t > 1);
if any(f)
	t(f) = 1;
	disp([' ## Out of bounds: ' int2str(length(f))])
end
t(1) = 0;
t(end) = 1;
t = t(:).';

if nargout > 0
	tick = t;
else
	if isdemo
		plot(x, pdf, '*', xi, yi, '-', xi, zi, '-', ...
				[t; t], [0*t; 0*t+1], 'g-')
		figure(gcf)
		s = [mfilename ' ' int2str(length(t)) ' ' int2str(length(pdf))];
		set(gcf, 'WindowButtonDownFcn', s)
		title(s)
	end
	assignin('caller', 'ans', t)
end
