function [x, y] = lambert_equal_area(lon, lat, doInverse, plon, plat, pang)

% lambert_equal_area -- Lambert equal-area projection.
%  [x, y] = lambert_eqarea(lon, lat, doInverse, plon, plat, pang)
%   returns the (x, y) positions corresponding to the Lambert
%   equal-area projection of (lon, lat), given in degrees.  The
%   mapping is centered on the given point (plon, plat, pang).
%   At the center point, the positive y-axis points upwards.
 
% svn $Id$
%=======================================================================
% Copyright (C) 1999 Dr. Charles R. Denham, ZYDECO.
%  All Rights Reserved.
%   Disclosure without explicit written consent from the
%    copyright owner does not constitute publication.
%=======================================================================
 
% Version of 27-Dec-1999 13:50:29.
% Updated    27-Dec-1999 13:50:29.

RCF = 180 / pi;

if nargin < 2, help(mfilename), return, end

if nargin < 3, doInverse = 0; end
if nargin < 4, plon = 0; end
if nargin < 5, plat = 0; end
if nargin < 6, pang = 0; end

if all(doInverse(:))
	xin = lon;
	yin = lat;
	[olon, olat] = inv_lambert_equal_area(xin, yin, plon, plat, pang);
	x = olon;
	y = olat;
	return
end

lon = lon / RCF;
lat = lat / RCF;

clat = cos(lat);
x = clat .* cos(lon);
y = clat .* sin(lon);
z = sin(lat);

if any([plat, plon, pang])
	[y, x] = rot1(y, x, plon);
	[z, x] = rot1(z, x, plat);
	[y, z] = rot1(y, z, pang);
end

h = sqrt(y.*y + z.*z);
ang = asin(h);
r = 2 * sin(ang / 2) / sqrt(2);
h(h == 0) = 1;
u = r .* y ./ h;
v = r .* z ./ h;

x = u/2;   % Normalize to 1 at 90 degrees radial distance.
y = v/2;

if any(pang)
	[x, y] = rot1(x, y, -pang);   % Restore North.
end


% ---------- inv_lambert_equal_area --------- %


function [lon, lat] = inv_lambert_equal_area(x, y, plon, plat, pang)

% inv_lambert_equal_area -- Inverse Lambert equal-area projection.

%  [lon, lat] = inv_lambert_equal_area(x, y, plon, plat, pang)
%   returns the (lon, lat) corresponding to the given Lambert
%   equal-area mapping (x, y), centered on the pole (plon, plat,
%   pang).
 
% Copyright (C) 1999 Dr. Charles R. Denham, ZYDECO.
%  All Rights Reserved.
%   Disclosure without explicit written consent from the
%    copyright owner does not constitute publication.
 
% Version of 27-Dec-1999 10:18:15.
% Updated    27-Dec-1999 10:18:15.

RCF = 180 / pi;

if any(pang)
	[x, y] = rot1(x, y, pang);
end

x = 2*x;
y = 2*y;
r = sqrt(x.*x + y.*y);
ang = 2*asin(r/2)*sqrt(2);

r(r == 0) = 1;

u = cos(ang);
v = sin(ang).*x./r;
w = sin(ang).*y./r;

x = u;
y = v;
z = w;

if any([plat, plon, pang])
	[y, z] = rot1(y, z, -pang);
	[z, x] = rot1(z, x, -plat);
	[y, x] = rot1(y, x, -plon);
end

lon = atan2(y, x) * RCF;
lat = asin(z) * RCF;


% ---------- rot1 --------- %


function [rx, ry] = rot1(x, y, deg)

% rot1 Planar rotation by an angle in degrees.
%  [rx, ry] = rot1(x, y, deg) rotates point X toward
%   Y by angle deg (in degrees).
%  ROT1 (no arguments) demonstrates itself.
 
% Copyright (C) 1992 Dr. Charles R. Denham, ZYDECO.
% All Rights Reserved.

% Version of 6-Jul-92 at 22:12:09.633.
% Updated    27-Dec-1999 10:16:20.

if nargin > 2
   xy = [x(:) y(:)].';
  else
   xy = x;
   deg = y;
end

rcf = 180 ./ pi;
rad = deg ./ rcf;
c = cos(rad); s = sin(rad);

r = [c -s; s c];

z = r * xy;

if nargout < 2
   rx = zeros(size(x));
   rx(:) = z;
else
   rx = zeros(size(x));
   ry = zeros(size(y));
   rx(:) = z(1, :); ry(:) = z(2, :);
end
