function zi = griddata1 (x, y, z, xi, yi, sr, ppq, expon)
					
% griddata1 -- 2-d data gridder.
%  griddata1(N, sr, ppq, expon) demonstrates itself.
%   Defaults: 100, inf, 1, and 1, respectively.
%  griddata1(x, y, z, xi, yi, sr, ppq, expon) is an inverse-
%   distance 2-d data gridder, which is guided by the search-
%   radius sr (default = inf), the points-per-quadrant ppq
%   (default = 1), and the exponent expon (default = 1).
%   The search-radius is the maximum distance permitted
%   between an interpolated position and the qualifying
%   data-points.  The points-per-quadrant is the maximum
%   number of (x, y) points to be used in each of the four
%   quadrants surrounding the corresponding (xi, yi) position.
%   The exponent for inverse-distance weighting defaults to 1,
%   denoting the simple reciprocal.  Invalid interpolates are
%   marked by NaN.
%  griddata1(x, y, z, xi, yi, [sr, ppq, expon]) is an
%   alternative syntax, in closer keeping with the "method"
%   syntax of the Matlab "griddata" function.

%  This began as John Evans' "gridgen_griddata" version
%   of a Rich Signell routine.  Cleaned-up by Denham to be
%   stand-alone, more versatile, and self-demonstrating.
 
% Modifications (Denham):
%
% svn $Id$
%=======================================================================
% Copyright (C) 2000 Dr. Charles R. Denham, ZYDECO.
%  All Rights Reserved.
%   Disclosure without explicit written consent from the
%    copyright owner does not constitute publication.
%=======================================================================
 
% Version of 20-Jul-2000 13:44:52.
% Updated    24-Jul-2000 11:12:17.

DEFAULT_SR = inf;    % Infinite search-radius.
DEFAULT_PPQ = 1;     % Nearest point in each quadrant.
DEFAULT_EXPON = 1;   % Very shallow unit-response.

if nargin < 1, x = 'demo'; help(mfilename), end

if isequal(x, 'demo'), x = 100; end

if isstr(x), x = eval(x); end

if length(x) == 1
	sr = DEFAULT_SR;
	ppq = DEFAULT_PPQ;
	expon = DEFAULT_EXPON;
	if nargin > 1, sr = y; end, if ischar(sr), sr = eval(sr); end
	if nargin > 2, ppq = z; end, if ischar(ppq), ppq = eval(ppq); end
	if nargin > 3, expon = xi; end
	if ischar(expon), expon = eval(expon); end
	n = x;
	x = rand(n, 1);
	y = rand(size(x));
	z = x+y;
	NGRID = 10;
	g = (0:NGRID)/NGRID;
	[xi, yi] = meshgrid(g, g);
	tic
	result = feval(mfilename, x, y, z, xi, yi, sr, ppq, expon);
	disp([' ## Elapsed time: ' int2str(round(toc)) ' s'])
	hold off
	try, set(gcf, 'Renderer', 'zbuffer'), catch, end
	surf(xi, yi, result)
	hold on
	plot3(x, y, z, '*k')
	hold off
	mfile = strrep(mfilename, '_', '\_');
	title([mfile ' ' int2str(n) ' ' num2str(sr) ' ' int2str(ppq) ...
			' ' num2str(expon)])
	xlabel x, ylabel y, zlabel z
	try, rotate3d on, catch, end
	return
end

if nargin  == 6 & length(sr) > 1
	for i = 1:length(sr), v{i} = sr(i); end
	zi = feval(mfilename, x, y, z, xi, yi, v{:});
	return
end

if nargin < 5, help(mfilename), return, end
if nargin < 6, sr = DEFAULT_SR; end
if nargin < 7, ppq = DEFAULT_PPQ; end
if nargin < 8, ppq = DEFAULT_EXPON; end

x = x(:); y = y(:); z = z(:);

search_radius = sr;
points_per_quadrant = max(ppq, 1);

n = length(xi(:));

% Loop through one grid-point at a time.

for k = 1:n
	
	% Locate the data lying within the search-radius,
	%  applied independently to x and y, rather than
	%  to actual radial distance.

	near_inds = find ( (abs(x-xi(k)) <= search_radius) & ...
		(abs(y-yi(k)) <= search_radius) );
	nearx = x(near_inds)-xi(k); 
	neary = y(near_inds)-yi(k);
	nearz = z(near_inds);
	
	quad1_inds = find ( (nearx>=0) & (neary>=0) );
	quad1x = nearx(quad1_inds); 
	quad1y = neary(quad1_inds);
	quad1z = nearz(quad1_inds);
	
	quad2_inds = find ( (nearx<0) & (neary>=0) );
	quad2x = nearx(quad2_inds); 
	quad2y = neary(quad2_inds);
	quad2z = nearz(quad2_inds);
	
	quad3_inds = find ( (nearx<0) & (neary<0) );
	quad3x = nearx(quad3_inds); 
	quad3y = neary(quad3_inds);
	quad3z = nearz(quad3_inds);
	
	quad4_inds = find ( (nearx>=0) & (neary<0) );
	quad4x = nearx(quad4_inds); 
	quad4y = neary(quad4_inds);
	quad4z = nearz(quad4_inds);
	
	% Select up to the maximum number per quadrant.

	if ( ~isempty(quad1_inds) )
		ranges_quad1 =  sqrt(quad1x.^2 + quad1y.^2);
		[sorted_range_quad1, range_inds] = sort(ranges_quad1);
		closest_quad1_inds = ...
			range_inds([1:1:min(length(range_inds),points_per_quadrant)]);
	else
		ranges_quad1 = [];
		closest_quad1_inds = [];
	end
	 
	if ( ~isempty(quad2_inds) )
		ranges_quad2 =  sqrt(quad2x.^2 + quad2y.^2);
		[sorted_range_quad2, range_inds] = sort(ranges_quad2);
		closest_quad2_inds = ...
			range_inds([1:1:min(length(range_inds),points_per_quadrant)]);
	else
		ranges_quad2 = [];
		closest_quad2_inds = [];
	end
	 
	if ( ~isempty(quad3_inds) )
		ranges_quad3 =  sqrt(quad3x.^2 + quad3y.^2);
		[sorted_range_quad3, range_inds] = sort(ranges_quad3);
		closest_quad3_inds = ...
			range_inds([1:1:min(length(range_inds),points_per_quadrant)]);
	else
		ranges_quad3 = [];
		closest_quad3_inds = [];
	end
	 
	if ( ~isempty(quad4_inds) )
		ranges_quad4 =  sqrt(quad4x.^2 + quad4y.^2);
		[sorted_range_quad4, range_inds] = sort(ranges_quad4);
		closest_quad4_inds = ...
			range_inds([1:1:min(length(range_inds),points_per_quadrant)]);
	else
		ranges_quad4 = [];
		closest_quad4_inds = [];
	end

	% Bundle all the qualifying data together.

	closest_ranges = [ranges_quad1(closest_quad1_inds); ...
					ranges_quad2(closest_quad2_inds); ...
					ranges_quad3(closest_quad3_inds); ...
					ranges_quad4(closest_quad4_inds); ];
				  
	closest_z = [quad1z(closest_quad1_inds); ...
					quad2z(closest_quad2_inds); ...
					quad3z(closest_quad3_inds); ...
					quad4z(closest_quad4_inds); ];

	% Use inverse-distance weighting with the nearest points
	%  in each quadrant. If none, the result is NaN.
	
	if ( isempty(closest_ranges) )
		zi(k) = NaN;
	else
		inv_ranges = closest_ranges .^ (-expon);
		zi(k) = sum ( ( inv_ranges ./ sum(inv_ranges) ) .* closest_z );
	end
		 
end

zi = reshape(zi, size(xi));
