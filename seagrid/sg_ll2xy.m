function [x, y] = sg_ll2xy(lon, lat)

% sg_ll2xy -- Stub for m_map-like "m_ll2xy" function.
%  [x, y] = sg_ll2xy(lon, lat) returns the (x, y) mappings
%   that correspond to geographic positions (lon, lat),
%   given in degrees.
 
% svn $Id$
%=======================================================================
% Copyright (C) 1999 Dr. Charles R. Denham, ZYDECO.
%  All Rights Reserved.
%   Disclosure without explicit written consent from the
%    copyright owner does not constitute publication.
%=======================================================================
 
% Version of 27-Dec-1999 13:37:03.
% Updated    25-Oct-2000 10:53:22.

[theProjection, theCenter] = sg_proj;
if isempty(theCenter)
	theCenter = [0 0 0];
end

plon = theCenter(1);
plat = theCenter(2);
pang = theCenter(3);

switch lower(theProjection)
case 'mercator'
	[x, y] = sg_mercator(lon, lat, 0, plon, plat, pang);
case 'transverse mercator'
	[x, y] = sg_mercator(lon, lat, 0, 0, 0, 90);
case 'stereographic'
	[x, y] = sg_stereographic(lon, lat, 0, plon, plat, pang);
case 'lambert equal area'
	[x, y] = sg_lambert_equal_area(lon, lat, 0, plon, plat, pang);
otherwise
	disp([' ## No such projection: ' theProjection])
end
