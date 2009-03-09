function [x, y] = mercator(lon, lat, doInverse, plon, plat, pang)

% mercator -- Mercator projection.
%  [x, y] = mercator(lon, lat, doInverse, plon, plat, pang)
%   returns the mercator projection of geographic positions
%   (lon, lat) in degrees, centered on the given point
%   (plon, plat, pang), whose default is (0, 0, 0).
%   If the doInverse flag is logically TRUE, the inverse
%   mapping is returned.  The Earth's radius is assumed
%   to be 1.  Regardless of the polar angle, positive y
%   on the map points toward North through the central
%   point.
%
%  For the transverse Mercator projection, use a pole
%  of (0, 0, 90).
 
% svn $Id$
%=======================================================================
% Copyright (C) 1999 Dr. Charles R. Denham, ZYDECO.
%  All Rights Reserved.
%   Disclosure without explicit written consent from the
%    copyright owner does not constitute publication.
%=======================================================================
 
% Version of 27-Dec-1999 10:18:15.
% Updated    27-Dec-1999 10:18:15.

RCF = 180 / pi;

if nargin < 2, help(mfilename), return, end

if nargin < 3, doInverse = 0; end
if nargin < 4, plon = 0; end
if nargin < 5, plat = 0; end
if nargin < 6, pang = 0; end

if all(doInverse(:))
	xin = lon;
	yin = lat;
	[olon, olat] = inv_mercator(xin, yin, plon, plat, pang);
	x = olon;
	y = olat;
	return
end

lon = lon / RCF;
lat = lat / RCF;

if any([plat, plon, pang])
	clat = cos(lat);
	x = clat .* cos(lon);
	y = clat .* sin(lon);
	z = sin(lat);
	[y, x] = rot1(y, x, plon);
	[z, x] = rot1(z, x, plat);
	[y, z] = rot1(y, z, pang);
	lon = atan2(y, x);
	lat = asin(z);
end

x = lon;
y = log(tan(pi/4 + lat/2));

if any(pang)
	[x, y] = rot1(x, y, -pang);   % Restore North.
end


% ---------- inv_mercator --------- %


function [lon, lat] = inv_mercator(x, y, plon, plat, pang)

% inv_mercator -- Inverse Mercator projection.

%  [lon, lat] = inv_mercator(x, y, plon, plat, pang) returns
%   the (lon, lat) corresponding to the given Mercator mapping
%   (x, y), centered on the pole (plon, plat, pang).
 
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

lon = x;
lat = (pi/2 - 2*atan(exp(-y)));

if any([plat, plon, pang])
	clat = cos(lat);
	x = clat .* cos(lon);
	y = clat .* sin(lon);
	z = sin(lat);
	[y, z] = rot1(y, z, -pang);
	[z, x] = rot1(z, x, -plat);
	[y, x] = rot1(y, x, -plon);
	lon = atan2(y, x);
	lat = asin(z);
end

lon = lon * RCF;
lat = lat * RCF;


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
