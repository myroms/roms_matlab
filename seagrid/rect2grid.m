function [result, error_norm] = rect2grid(z, zrect, theCorners, theSize)

% rect2grid -- Orthogonal grid from RECT result.
%  [result, error_norm] = rect2grid(z, zrect, theCorners, theSize)
%   produces a curvilinear orthogonal grid by interpolating the
%   complex contour z, using zrect, the result of applying the
%   conformal mapper RECT to z for theCorners (indices).  If zrect
%   is empty, the RECT routine is called to compute it.  If zrect
%   is a scalar, that number of RECT iterations will be performed.
%   The size of the w output grid (complex matrix), including the
%   perimeter, is determined by theSize.  The returned error_norm
%   is the norm of the respective laplacians.
%  rect2grid(nPoints) demonstrates itself with a random z contour
%   of nPoints (default = 20).
 
% svn $Id$
%=======================================================================
% Copyright (C) 1998 Dr. Charles R. Denham, ZYDECO.
%  All Rights Reserved.
%   Disclosure without explicit written consent from the
%    copyright owner does not constitute publication.
%=======================================================================
 
% Version of 21-Oct-1998 20:50:16.
% Updated    23-Jun-2000 14:29:17.

% Note: The splined interpolations tend to be too wild
%  for general use here.  The safest method is to sample
%  the z curve at very fine intervals, then use linear
%  interpolation.
%
% Note: Presently, we apply rect to the existing boundary
%  at equal-intervals along the edges, then interpolate
%  to even-intervals in the rect domain, then solve the
%  Poisson system, then interpolate back to the original
%  spacing.  It might make better sense to iteratively
%  apply rect, then resample the original boundary, until
%  the rect-domain itself shows more or less equal-spacing
%  along all edges.  We would adjust one edge before moving
%  to the next.  Then, we would solve the Poisson system,
%  then interpolate backwards as before.  The aim is to avoid
%  crowding of mapped points in the rect-domain, where things
%  are often not very smooth.  If we have to interpolate,
%  let it be restricted to the original domain.
		
% Note: The rectangular mapping itself is based on a boundary
%  composed of linear segments, which implies that any subsequent
%  interpolations within the Poisson solution-grid should be done
%  bi-linearly, rather than via splines.  (The bilinear interpolants
%  of a rectangular grid have Laplacian = 0 automatically.)  By
%  extension, the best way to interpolate a Poisson grid is
%  linearly if we are taking the FFT approach.

% Note: We should think about an iterative scheme is which the
%  points to be mapped are adjusted until the mapping yields
%  a more-or-less equally-spaced distribution.  It will be
%  tricky, because interated-splines are badly behaved.

if (0)
	theInterpFcn = 'splineq';   % Interpolation function.
	theInterpMethod = 1;   % interpolation method.
else
	theInterpFcn = 'interp1';   % Interpolation function.
	theInterpMethod = 'linear';   % interpolation method.
end

if nargin < 1, help(mfilename), z = 'demo'; end
if isequal(z, 'demo'), z = 20; end
if ischar(z), z = eval(z); end

if length(z) == 1
	n = max(z, 4);
	z = rand(n, 1) + sqrt(-1)*rand(n, 1);
	z = z - mean(z);
	a = angle(z);
	[a, i] = sort(a);
	jitter = 0.1;
	z = (1 + jitter * rand(size(a))) .* exp(sqrt(-1)*a);
	[ignore, nn] = sort(rand(1, length(z)-1));
	theCorners = sort([1 nn(1:3)+1]);
	
	k = 1; ktest = n*sqrt(n);
	while k < ktest
		k = 2*k;
	end
	theSize = [k k] + 1;
	
	ztemp = z;
	ztemp(end+1) = z(1);
	z = [];
	c = [theCorners length(ztemp)];
	theCorners = [];
	
	zplot = [];
	
	for i = 1:4
		j = c(i):c(i+1);
		[ji, pp] = splineq(j, ztemp(j), n, 1);
		zitemp = ppval(pp, ji);
		zitemp(end) = [];
		theCorners = [theCorners length(z)+1];
		z = [z zitemp];
	end
	
	if (0)   % Set to 1 to see orthogonality.
		z = rect(z, ceil(sqrt(n)), theCorners);
	end
	
	tic

	FLAG = 0;
	
	zrect = [];
	
	if FLAG
		zrect = rect(z, ceil(sqrt(n)), theCorners);
	end
	
	[w, err] = rect2grid(z, zrect, theCorners, ceil(theSize/2));
	
	disp([' ## Elapsed time: ' num2str(toc) ' seconds.'])
	
	if ~isempty(w)
		u = real(w); v = imag(w);
		u_err = real(err); v_err = imag(err);
		subplot(1, 1, 1)
		if FLAG, subplot(1, 2, 1), end
		hold off
		h = plot(u, v, 'g-', u.', v.', 'b-');
		hold on
		z(end+1) = z(1);
		x = real(z); y = imag(z);
		plot(x, y, 'bo', ...
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
		if FLAG
			subplot(1, 2, 2)
			plot(real(zrect), imag(zrect), '-o')
			axis equal
			zoomsafe 0.9, zoomsafe
		end
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
			
hasMex = 0;   % Note this override.

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
	tolerance = 0.00001;
	
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
		
% Verbose quality of mapping.

	if (0)
		hello rect2grid line 150
		non_straight_percent = 100*(1-straight)
		d13 = abs(diff(zrect(theCorners([1 3]))));
		d24 = abs(diff(zrect(theCorners([2 4]))));
		non_rect_percent = 100*2*abs((d13-d24)/(d13+d24))
		fig = gcf;
		f = findobj('Type', 'figure', 'Tag', mfilename);
		if isempty(f)
			f = figure('Name', mfilename', 'Tag', mfilename);
		end
		figure(f)
		plot(real(zrect), imag(zrect), '-o', ...
			real(zrect(theCorners(1))), imag(zrect(theCorners(1))), '*')
		set(gca, 'XLim', [-0.1 1.1]*max(real(zrect)), ...
				'YLim', [-0.1 1.1]*max(imag(zrect)))
		figure(fig)
		drawnow
	end

	if norm(1-straight) > tolerance
		disp([' ## rect2grid: Conformal mapping not successful'])
		disp(['               after ' int2str(ntimes) ' iterations.'])
		if nargout > 0, result = []; error_norm = []; end
		return
	end
end

% Make zrect perfectly rectangular: not necessary.

if (0)
	x = real(zrect); x = [x(:); x(1)];
	y = imag(zrect); y = [y(:); y(1)];
	c = [theCorners(:); length(x)].';
	y(c(1):c(2)) = 0;
	x(c(2):c(3)) = 1;
	y(c(3):c(4)) = mean(y(c(3):c(4)));
	x(c(4):c(5)) = 0;
	zrect = x(1:end-1) + sqrt(-1)*y(1:end-1);
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

rdist = [0; cumsum(abs(diff([zrect(:); zrect(1)])))];
rdist = rdist - min(rdist); rdist = rdist / max(rdist);

z(end+1) = z(1);
x = real(z);
y = imag(z);

c = [theCorners(:); length(x)].';   % Physical and mapped corners.
d = cumsum([1 m-1 n-1 m-1 n-1]);    % Corners around the matrix.

xi = zeros(size(ind));
yi = zeros(size(ind));

for i = 1:4
	j = c(i):c(i+1);
	k = d(i):d(i+1);
	rd = rdist(j); rd = rd - min(rd); rd = rd / max(rd);
	xx = x(j);
	yy = y(j);
	dd = k; dd = dd - min(dd); dd = dd / max(dd);
	
	TESTING = 1;   % <== NOTE.
	TESTING = 0;   % <== NOTE.
	
	if TESTING
		[ki, pp] = splineq(j, z(j), length(k), 1);
		zi = ppval(pp, ki);
		xi(k) = real(zi);
		yi(k) = imag(zi);
	else
		xi(k) = feval(theInterpFcn, rd, xx, dd, theInterpMethod).';
		yi(k) = feval(theInterpFcn, rd, yy, dd, theInterpMethod).';
	end
	
end

if (0)
	hello(mfilename)
	c, y, d = diff(y);
end

% Sprinkle interpolated values along the perimeter.

u = zeros(theSize); v = zeros(theSize);
u(ind) = xi; v(ind) = yi;

% Aspect ratio of the rectangle.

if (1)
	dx = 1; dy = 1;   % Square.
else
	dx = abs(zrect(2)-zrect(1)) / m;   % Rectangle.
	dy = abs(zrect(3)-zrect(2)) / n;   % Rectangle.
end

% Solve Laplace's equation inside the boundary.  At this
%  stage, it is slightly advantageous to use arrays u and
%  v whose sizes are a power-of-two plus one.

if ~hasMex
	isSlope = 0;
	u = fps(u, isSlope, dx, dy);
	v = fps(v, isSlope, dx, dy);
else   % Use MEXSEPELI.
	l2 = m-1;
	m2 = n-1;
	seta = linspace(0, 1, n);
	sxi = linspace(0, 1, m);
	[u, v] = feval('mexsepeli', u, v, l2, m2, seta, sxi);
end

w = u + sqrt(-1).*v;

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
