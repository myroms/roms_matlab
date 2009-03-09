function [theProj, theCenter] = sg_proj(theProjection, theProjectionCenter)

% sg_proj -- Stub for m_map-like "m_proj" function.
%  sg_proj(theProjection) sets and/or gets the current
%   map projection: {'Mercator', 'Transverse Mercator', 'Stereographic'}.
%  [theProj, theCenter] = sg_proj(theProjection, theProjectionCenter)
%   also sets/gets the center-of-projection (lon, lat, angle), given
%   in degrees.
 
% svn $Id$
%=======================================================================
% Copyright (C) 1999 Dr. Charles R. Denham, ZYDECO.
%  All Rights Reserved.
%   Disclosure without explicit written consent from the
%    copyright owner does not constitute publication.
%=======================================================================
 
% Version of 27-Dec-1999 11:26:45.
% Updated    27-Dec-1999 13:37:03.

persistent CURRENT_PROJECTION
persistent CURRENT_PROJECTION_CENTER

if nargin > 0
	switch lower(theProjection)
	case {'mercator', 'stereographic', 'lambert equal area'}
		CURRENT_PROJECTION = theProjection;
	otherwise
		disp([' ## No such projection: ' theProjection])
	end
end

if nargin > 1
	CURRENT_PROJECTION_CENTER = theProjectionCenter;
end

if nargout > 0
	theProj = CURRENT_PROJECTION;
	theCenter = CURRENT_PROJECTION_CENTER;
end
