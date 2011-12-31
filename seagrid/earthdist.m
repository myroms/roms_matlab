function theResult = earthdist(alon, alat, blon, blat, radius)

% earthdist -- Distance in meters between two lon/lats.
%  earthdist(alon, aloat, blon, blat, radius) returns the
%   distance in maters between locations (alon, alat)
%   and (blon, blat).  The default earth radius is
%   assumed to be 6371.315*1000 meters, the radius for
%   a sphere of equal-volume.

% svn $Id$
%=======================================================================
% Copyright (C) 1996 Dr. Charles R. Denham, ZYDECO.
% All Rights Reserved.
%   Disclosure without explicit written consent from the
%    copyright owner does not constitute publication.
%=======================================================================

%  This function was modified by Hernan G. Arango (11/25/11) so the
%  great circle distance is identical to the formula used in grid
%  refinement.

if nargin < 4, help(mfilename), return, end
if nargin < 5, radius = 6371.315*1000; end   % meters.

RCF = 180 / pi;

alon = alon / RCF;
alat = alat / RCF;
blon = blon / RCF;
blat = blat / RCF;

alpha = sin(alat).*sin(blat) + cos(alat).*cos(blat).*cos(blon-alon);

ind = find (abs(alpha)>1);
if (~isempty(ind)), alpha(ind) = sign(alpha(ind)); end;

alpha = acos(alpha);
result = radius .* alpha;

if nargout > 0
	theResult = result;
else
	disp(result)
end
