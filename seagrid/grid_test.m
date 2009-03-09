function grid_test(m, n)

% grid_test -- Grid generation demo.
%  grid_test(m, n) constructs a grid of size [m n].  This
%   is a demonstration routine which attempts to lay out
%   the algorithm in an intuitive manner.
 
% svn $Id$
%=======================================================================
% Copyright (C) 1999 Dr. Charles R. Denham, ZYDECO.
%  All Rights Reserved.
%   Disclosure without explicit written consent from the
%    copyright owner does not constitute publication.
%=======================================================================
 
% Version of 27-May-1999 22:22:17.
% Updated    03-Jun-1999 15:09:18.

if nargin < 1, m = 11; end
if nargin < 2, n = m; end

if ischar(m), m = eval(m); end
if ischar(n), n = eval(n); end

x = [1 4 3 1];
y = [1 1 2 2];

z = x + sqrt(-1)*y;
c = [1 2 3 4];   % corner indices.

z(end+1) = z(1);
c(end+1) = length(z);

counts = [m n m n];

zz = [];
cc = [];
for i = 1:4
	j = c(i):c(i+1);
	[ji, pp] = splineq(j, z(j), counts(i), 1);
	zi = ppval(pp, ji);
	zi(end) = [];
	cc = [cc length(zz)+1];
	zz = [zz zi(:).'];
end

z = zz; c = cc;

subplot(2, 2, 1)
x = real(z); y = imag(z);
plot(x, y, '+-', x(1), y(1), 'o')
title z, xlabel x, ylabel y

zoomsafe 0.9, zoomsafe

hasMex = (exist('mexrect', 'file') == 3);

iteration = 0;
straightness = 0;
zrect = z;
while abs(1-straightness) > 0.0001 & iteration < max(m, n)
	if ~hasMex
		zrect = rect(zrect, 1, c(1), c(2), c(3), c(4));
	else
		zrect = mexrect(zrect, length(zrect), c(1), c(2), c(3), c(4));
	end
	ztemp = zrect;
	ztemp(end+1) = ztemp(1);
	ctemp = c;
	ctemp(end+1) = length(ztemp);
	d = sum(abs(diff(ztemp)));
	dc = sum(abs(diff(ztemp(ctemp))));
	straightness = dc ./ d;
	iteration = iteration+1;
	disp([' ## iteration ' int2str(iteration) ...
		'; straightness = ' num2str(100*straightness) ' percent'])
end

subplot(2, 2, 2)
xr = real(zrect); yr = imag(zrect);
plot(xr, yr, '+-', xr(1), yr(1), 'o')
title zrect, xlabel u, ylabel v
zoomsafe 0.9, zoomsafe

figure(gcf)
set(gcf, 'Name', mfilename)

ind = zeros(m, n);
ind(:) = 1:prod(size(ind));
index = [ind(1:end-1, 1).', ind(end, 1:end-1), ...
			ind(end:-1:2, end).', ind(1, end:-1:2)];

z = z(:).';
z(end+1) = z(1);

zrect = zrect(:).';
zrect(end+1) = zrect(1);
c(end+1) = length(z);

d = [0 cumsum(abs(diff(z)))];		% Original boundary.
e = [0 cumsum(abs(diff(zrect)))];	% Mapped boundary.

d12 = d(c(1):c(2));
d12 = d12 - min(d12); d12 = d12 ./ max(d12);
d23 = d(c(2):c(3));
d23 = d23 - min(d23); d23 = d23 ./ max(d23);

e12 = e(c(1):c(2));
e12 = e12 - min(e12); e12 = e12 ./ max(e12);
e23 = e(c(2):c(3))
e23_scale = max(e23) - min(e23)
e23 = e23 - min(e23); e23 = e23 ./ max(e23);

hasMex = (exist('mexsepeli', 'file') == 3)

theSlopeFlag = 1;
theSlopeFlag = 0;

zz = zeros(size(z));

if ~hasMex
	for i = 1:4
		j = c(i):c(i+1);
		dd = d(j);
		dd = dd - min(dd); dd = dd ./ max(dd);
		ee = e(j);
		ee = ee - min(ee); ee = ee ./ max(ee);
		pp = splinep(ee, z(j), theSlopeFlag);
		if i == 1
			dd = d12;
		elseif i == 2
			dd = d23;
		elseif i == 3
			dd = fliplr(1-d12);
		elseif i == 4
			dd = fliplr(1-d23);
		end
		zz(j) = ppval(pp, dd);
	end
	zz(end) = [];
	u = zeros(m, n);
	v = zeros(m, n);
	u(index) = real(zz);
	v(index) = imag(zz);
	u = fps(u);
	v = fps(v);
else   % Use Mex-file.
	seta = e12;
	sxi = e23 * e23_scale;   % Note.
	for i = 1:4
		j = c(i):c(i+1);
		ee = e(j);
		ee = ee - min(ee); ee = ee ./ max(ee);
		if i == 1
			pp = splinep(ee, z(j), theSlopeFlag);
			zz(j) = ppval(pp, e12);
		elseif i == 2
			pp = splinep(ee, z(j), theSlopeFlag);
			zz(j) = ppval(pp, e23);
		elseif i == 3
			pp = splinep(ee, z(j), theSlopeFlag);
			zz(j) = ppval(pp, fliplr(1-e12));
		elseif i == 4
			pp = splinep(ee, z(j), theSlopeFlag);
			zz(j) = ppval(pp, fliplr(1-e23));
		end
	end
	zz(end) = [];
	u = zeros(m, n);
	v = zeros(m, n);
	u(index) = real(zz);
	v(index) = imag(zz);
	[u, v] = mexsepeli(u, v, m-1, n-1, sxi, seta);   % Note calling sequence.
end

uu = del2(u);
vv = del2(v);
laplacian = norm(uu(2:end-1, 2:end-1)) + norm(vv(2:end-1, 2:end-1))

subplot(2, 1, 2)
plot(u, v, u.', v.')
title gridded, xlabel x, ylabel y
axis equal
zoomsafe 0.9, zoomsafe
