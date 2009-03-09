function [result, error_norm] = mexrect2grid(z, zrect, theCorners, theSize)

% mexrect2grid -- Orthogonal grid from RECT result via mex-files.
%  [result, error_norm] = mexrect2grid(z, zrect, theCorners, theSize)
%   produces a curvilinear orthogonal grid by interpolating the
%   complex contour z, using zrect, the result of applying the
%   conformal mapper RECT to z for theCorners (indices).  If zrect
%   is empty, the RECT routine is called to compute it.  If zrect
%   is a scalar, that number of RECT iterations will be performed.
%   The size of the w output grid (complex matrix), including the
%   perimeter, is determined by theSize.  The returned error_norm
%   is the norm of the respective laplacians.  This routine uses
%   the "mexrect" and "mexsepeli" mex-files if available; otherwise,
%   it calls the "rect" and "fps" m-files.
%  mexrect2grid(nPoints) demonstrates itself with a random z contour
%   of nPoints (default = 20).
 
% svn $Id$
%=======================================================================
% Copyright (C) 1998 Dr. Charles R. Denham, ZYDECO.
%  All Rights Reserved.
%   Disclosure without explicit written consent from the
%    copyright owner does not constitute publication.
%=======================================================================
 
% Version of 21-Oct-1998 20:50:16.
% Updated    09-Jun-1999 10:24:25.

if nargin < 1, help(mfilename), z = 'demo'; end
if isequal(z, 'demo'), z = 20; end
if ischar(z), z = eval(z); end

if length(z) == 1
	n = z;
	z = rand(n, 1) + sqrt(-1)*rand(n, 1);
	z = z - mean(z);
	a = angle(z);
	[a, i] = sort(a);
	jitter = 0.1;
	z = (1 + jitter * rand(size(a))) .* exp(sqrt(-1)*a);
	[ignore, nn] = sort(rand(1, length(z)-1));
	theCorners = sort([1 nn(1:3)+1]);
	theSize = 2*[n n];
	zrect = [];
	tic
	[w, err] = feval(mfilename, z, zrect, theCorners, ceil(theSize/2));
	disp([' ## Elapsed time: ' num2str(toc) ' seconds.'])
	if ~isempty(w)
		u = real(w); v = imag(w);
		u_err = real(err); v_err = imag(err);
		x = real(z); y = imag(z);
		x = [x; x(1)]; y = [y; y(1)];  % Make closed curve.
		u1 = u(:, 2:end-1);   % Trim the grid.
		v1 = v(:, 2:end-1);
		u2 = u(2:end-1, :).';
		v2 = v(2:end-1, :).';
		h = plot(u1, v1, 'g-', u2, v2, 'b-');
		hold on
		plot(x, y, 'r-', ...
				x(theCorners), y(theCorners), 'ro', ...
				x(theCorners(1)), y(theCorners(1)), 'r*')
		hold off
		xlabel('x'), ylabel('y')
		theCommand = [mfilename ' ( ' int2str(n) ' )'];
		title(theCommand)
		set(gcf, 'ButtonDownFcn', theCommand)
		figure(gcf)
		axis equal
		zoomsafe 0.9, zoomsafe
	end
	error_norm = [real(err) imag(err)];
	if nargout > 0
		result = w;
	else
		disp([' ## error_norm = ' sprintf('%0.4g %0.4gi', error_norm)])
	end
	return
end

% Check for mex-files "mexrect" and "mexsepeli".

hasMex = (exist('mexrect', 'file') == 3) & ...
			(exist('mexsepeli', 'file') == 3);

% If no actual "zrect" is given, apply RECT until
%  the straightness of the result deviates from
%  1.0 by no more than 0.1 percent.

if length(zrect) < 2
	if length(zrect) == 1
		ntimes = zrect;
	else
		ntimes = ceil(sqrt(length(z)));
	end
	zrect = z(:);
	tolerance = 0.001;
	for i = 1:ntimes
		if ~hasMex
			[zrect, straight] = feval('rect', zrect, 1, theCorners);
		else
			zrect = feval('mexrect', zrect, length(zrect), ...
				theCorners(1), theCorners(2), theCorners(3), theCorners(4));
			ztemp = zrect;
			ztemp(end+1) = ztemp(1);
			ctemp = theCorners;
			ctemp(end+1) = ctemp(1);
			num = sum(abs(diff(ztemp(ctemp))));
			den = sum(abs(diff(ztemp)));
			straight = num./den;
		end
		if norm(1-straight) <= tolerance, break, end
		disp([' ## RECT Iteration ' int2str(i) ...
				': straightness = ' num2str(straight*100) ' percent.'])
	end
	if norm(1-straight) > tolerance
		disp([' ## rect2grid: Conformal mapping not successful'])
		disp(['               after ' int2str(ntimes) ' iterations.'])
		if nargout > 0, result = []; error_norm = []; end
		return
	end
end

% Desired size.

if length(theSize) == 1, theSize = theSize * [1 1]; end

m = theSize(1); n = theSize(2);

% Get indices of matrix perimeter.

temp = zeros(theSize);
temp(:) = 1:prod(theSize);

ind = [];
ind = [ind; temp(1:m-1, 1)];
ind = [ind; temp(m, 1:n-1).'];
ind = [ind; temp(m:-1:2, n)];
ind = [ind; temp(1, n:-1:1).'];

% Interpolate around the "zrect" boundary
%  as a function of distance along the physical
%  boundary.

zrect = zrect(:).';
zrect(end+1) = zrect(1);

rdist = [0 cumsum(abs(diff(zrect)))];
rdist = rdist - min(rdist); rdist = rdist / max(rdist);

z = z(:).';
z(end+1) = z(1);
c = theCorners;
c(end+1) = length(z);

d = cumsum([1 m-1 n-1 m-1 n-1]);    % Corners around the matrix.

zi = zeros(size(ind));

slopeFlag = 1;

for i = 1:4
	j = c(i):c(i+1);   % Data corners.
	k = d(i):d(i+1);   % Matrix corners.
	rd = rdist(j); rd = rd - min(rd); rd = rd / max(rd);
	if i == 1
		pp = splinep(linspace(0, 1, length(rd)), rd, slopeFlag);
		rd1 = ppval(pp, linspace(0, 1, m));
		pp = splinep(rd, z(j), slopeFlag);
		zi(k) = ppval(pp, rd1);
	elseif i == 2
		pp = splinep(linspace(0, 1, length(rd)), rd, slopeFlag);
		rd2 = ppval(pp, linspace(0, 1, n));
		pp = splinep(rd, z(j), slopeFlag);
		zi(k) = ppval(pp, rd2);
	elseif i == 3
		pp = splinep(rd, z(j), slopeFlag);
		zi(k) = ppval(pp, fliplr(1 - rd1));
	elseif i == 4
		pp = splinep(rd, z(j), slopeFlag);
		zi(k) = ppval(pp, fliplr(1 - rd2));
	end
end

% Sprinkle interpolated values along the perimeter.

u = zeros(theSize); v = zeros(theSize);
u(ind) = real(zi); v(ind) = imag(zi);

% Aspect ratio of the rectangle.

if (1)
	dx = 1; dy = 1;   % Square.
else
	dx = abs(zrect(2)-zrect(1)) / m;   % Rectangle.
	dy = abs(zrect(3)-zrect(2)) / n;   % Rectangle.
end

% Solve Laplace's equation inside the boundary.

if ~hasMex
	isSlope = 0;
	u = feval('fps', u, isSlope, dx, dy);
	v = feval('fps', v, isSlope, dx, dy);
else   % Use MEXSEPELI.
	l2 = m-1;
	m2 = n-1;
	seta = rdist(c(2):c(3));
	sxi = rdist(c(1):c(2));
	[u, v] = feval('mexsepeli', u, v, l2, m2, seta, sxi);
end

w = u + sqrt(-1)*v;

if nargout > 0
	result = w;
else
	disp(w)
end

if nargout > 1
	del2_u = 4*del2(u);
	err_norm_u = norm(del2_u(2:end-1, 2:end-1));
	del2_v = 4*del2(v);
	err_norm_v = norm(del2_v(2:end-1, 2:end-1));
	error_norm = err_norm_u + sqrt(-1).*err_norm_v;
end
