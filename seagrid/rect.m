function [theResult, theStraightness] = rect(z, ntimes, n1, n2, n3, n4)

% rect -- Map contour points onto a rectangle.
%  rect(z, ntimes, n1, n2, n3, n4) maps the complex z
%   contour points, in counter-clockwise sequence,
%   onto a unit-rectangle, using ntimes iterations.
%   The n1..n4 are the indices of the four
%   corner-points.
%  [theResult, theStraightness] = rect(...) also
%   returns a measure of straightness, which
%   should be very close to 1.  Convergence can be
%   checked by computing the norm of the difference
%   of two successive iterations.
%  rect(n, ntimes) demonstrates itself with n random
%   points (default is 40) and ntimes iterations
%   (default is sqrt(n)).  The corner-positions
%   are chosen randomly.
%
% Reference: Ives, D.D. and R.M. Zacharias, Conformal
%  mapping and orthogonal grid generation, Paper No.
%  87-2057, AIAA/SAE/ASME/ASEE 23rd Joint Propultion
%  Conference, San Diego, California, June, 1987.

% The present routine is modified slightly to avoid
%  zero-divides.
 
% svn $Id$
%=======================================================================
% Copyright (C) 1997 Dr. Charles R. Denham, ZYDECO.
%  All Rights Reserved.
%   Disclosure without explicit written consent from the
%    copyright owner does not constitute publication.
%=======================================================================
 
% Version of 05-Jun-1998 09:07:06.
% Updated    20-Oct-1998 22:21:56.
% Revised    26-Oct-1998 14:29:11.

if nargin < 1, help(mfilename), z = 'demo'; end

if isequal(z, 'demo'), z = 40; end

if isstr(z), z = eval(z); end

if length(z) == 1
	if nargin < 2, ntimes = ceil(sqrt(z)); end
	if ischar(ntimes), ntimes = eval(ntimes); end
	n = fix((z+3)/4) .* 4;
	x = rand(1, n) - 0.5;
	y = rand(1, n) - 0.5;
	z = x + sqrt(-1) .* y;
	[a, ind] = sort(angle(z));
	z = z(ind);
	z = z + exp(sqrt(-1).*a);
	z = 10 * z;
	if (0)
		nn = round([1 n/4 2*n/4 3*n/4]);
	else   % Random corner-points.
		[s, index] = sort(rand(1, length(z)-1));
		nn = [1 (sort(index(1:3))+1)];
	end
	delete(get(gcf, 'Children'))
	subplot(1, 2, 1)
	plot(real(z), imag(z), '-+', real(z(nn)), imag(z(nn)), 'ro')
	title('Original')
	axis equal
	drawnow
	figure(gcf)
	et = 0;
	result = z;
	for iter = 1:ntimes
		z = result;
		tic
		[result, straightness] = rect(z, 1, nn);
		if ~isfinite(straightness), break, end
		et = et + toc;
		subplot(1, 2, 2)
		plot(real(result), imag(result), '-+', ...
				real(result(nn)), imag(result(nn)), 'ro')
		title(['Mapped by RECT: iteration #' int2str(iter)])
		axis equal
		figure(gcf)
		drawnow
	end
	title(['Mapped by RECT'])
	disp([' ## Elapsed time: ' int2str(et) ' seconds.'])
	if nargout > 0
		theResult = result;
		theStraightness = straightness;
	else
		disp([' ## Straightness: ' num2str(straightness)])
	end
	zoomsafe
	s = [mfilename ' ' int2str(n) ' ' int2str(ntimes)];
	set(gcf, 'ButtonDownFcn', s)
	return
end

z = z(:);
if z(1) == z(end), z(end) = []; end
n = length(z);

if nargin < 2, ntimes = 1; end
if nargin == 3 & length(n1) == 4
	n4 = n1(4);
	n3 = n1(3);
	n2 = n1(2);
	n1 = n1(1);
else
	if nargin < 3, n1 = 1; end
	if nargin < 4, n2 = n1 + fix(n/4); end
	if nargin < 5, n3 = n2 + fix(n/4); end
	if nargin < 6, n4 = n3 + fix(n/4); end
end

nn = [n1 n2 n3 n4];

r = zeros(size(z));
t = zeros(size(z));

track_errors = 0;

for iter = 1:ntimes
	
	if ntimes > 10
		disp([' ## Iterations remaining: ' int2str(ntimes-iter+1)])
	end

	z_old = z;
	
	for i = 1:n   % The big loop.
		if rem(n-i+1, 100) == 0 & 0
			disp([' ## Remaining: ' int2str(n-i+1)])
		end
		im = n - rem(n-i+1, n);	% Index of previous point.
		ip = 1 + rem(i, n);		% Index of next point.
		zd = 1;
		if z(ip) ~= z(i)
			zd = (z(im) - z(i)) ./ (z(ip) - z(i));
		else
			zd = (z(im) - z(i)) ./ sqrt(eps);
		end
		alpha = atan2(imag(zd), real(zd));   % Current angle.
		if alpha < 0
			alpha = alpha + 2 .* pi;
		elseif alpha == 0
			alpha = sqrt(eps);
		end
		pwr = pi ./ alpha;   % Power for pi outcome.
		if any(i == nn)   % Corners.
			if z(im) ~= z(i)
				zd = (z(im) - z(i)) ./ abs(z(im) - z(i));   % Unit-vector.
			else
				zd = 1;
			end
			z = sqrt(-1) .* z ./ zd;
			pwr = pwr ./ 2;   % Power for pi/2 outcome.
		end
	
% Compute the unwrapped phases.
		
		j = 2:n;
		zd = z(rem(j+i-2, n)+1) - z(i);
		r(j) = abs(zd);
		t(j) = atan2(imag(zd), real(zd)) - 6.*pi;
		if (1)
			t = unwrap1(t);   % Local subroutine.
		else
			t = unwrap(t);   % Matlab subroutine.
		end
		temp = t(2:n).*pwr;
		
		pmax = max(temp);
		pmin = min(temp);
		if pmax ~= pmin
			pwr = min(pwr, 1.98 .* pi .* pwr ./ (pmax - pmin));
		end
		z(i) = 0;
		
		z(rem(j+i-2, n)+1) = r(j).^pwr.*exp(sqrt(-1).*t(j).*pwr);
		zd = 1 ./ (z(n2) - z(n1));  % Possible zero-divide.
		z0 = z(n1);
		z = (z - z0) .* zd;
	end

% Compute the straightness of the perimeter
%  by comparing the perimeter along the data
%  with the perimeter just around the corners.
%  The mapping may appear to be straight
%  without actually being rectangular.

	p1 = sum(abs(diff([z(nn); z(nn(1))])));
	p2 = sum(abs(diff([z; z(1)])));
	straightness = p1 ./ p2;

% Norm of changes as measure of convergence.
%  The mapping may converge without actually
%  producing a rectangular shape.

	change = 2 .* norm(z - z_old) ./ norm(z + z_old);

% Progress report.

	s = round(1000 .* straightness) ./ 10;
	c = round(1000 .* change) ./ 10;
	if nargout < 2 & 0
		disp([' ## ' num2str(s) ' % straight; ' num2str(c) ' % change.'])
	end
	
end

if track_errors > 0   % Not ever checked.
	disp([' ## Tracking errors: ' int2str(track_errors)])
end

if nargout > 0
	theResult = z;
	theStraightness = straightness;
else
	disp(z), straightness
end


function theResult = unwrap1(p)

% unwrap1 -- Unwrap radian phases modulo 2*pi.
%  unwrap1(p) performs a modulo-2*pi unwrapping
%   of the phases of radian values p, working
%   columnwise if p is a matrix.
%  unwrap1 (no argument) demonstrates itself.
 
% Copyright (C) 1998 Dr. Charles R. Denham, ZYDECO.
%  All Rights Reserved.
%   Disclosure without explicit written consent from the
%    copyright owner does not constitute publication.
 
% Version of 09-Nov-1998 11:02:40.

if nargin < 1, help(mfilename), p = 'demo'; end

if isequal(p, 'demo'), p = 24; end

if ischar(p), p = eval(p); end

if length(p) == 1
	m = p;
	p = rand(m, 1)*2*pi;
	r = unwrap1(p);
	t = (1:m).';
	plot(t, p, ':o', t, r, '-o')
	legend('original', 'unwrapped')
	title('unwrap1 demo')
	xlabel('time'), ylabel('phase (radians)')
	figure(gcf)
	return
end

[oldM, n] = size(p);
if oldM == 1, p = p.'; end
[m, n] = size(p);

mn = ones(m, 1) * min(p);   % Column-wise minima.
p = rem(p-mn, 2*pi) + mn;   % Flatten to 2*pi range.
d = [p(1, :); diff(p)];   % Differences.
up = (d > pi); down = (d <= -pi);   % Find big jumps.
c = cumsum(down-up)*2*pi;   % Corrections.
result = p+c;

if oldM == 1, result = result.'; end

if nargout > 0
	theResult = result;
else
	disp(result)
end
