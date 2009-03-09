function [lon, lat] = sg_xy2ll(x, y)

% sg_ll2xy -- Stub for m_map-like "m_ll2xy" function.
%  [lon, lat] = sg_xy2ll(x, y) returns the (lon, lat)
%   geographic positions, in degrees, that correspond
%   to the given (x, y).
 
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
	[lon, lat] = sg_mercator(x, y, 1, plon, plat, pang);
case 'transverse mercator'
	[lon, lat] = sg_mercator(x, y, 1, 0, 0, 90);
case 'stereographic'
	[lon, lat] = sg_stereographic(x, y, 1, plon, plat, pang);
case 'lambert equal area'
	[lon, lat] = sg_lambert_equal_area(x, y, 1, plon, plat, pang);
otherwise
	disp([' ## No such projection: ' theProjection])
end
