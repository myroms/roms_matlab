function help_projection(self)

%    SeaGrid offers "Mercator" (default) and "Stereographic"
% conformal mapings, of which the Mercator projection looks
% the best here.  The SeaGrid "Setup" dialog allows the
% projection, center-of-projection, and geographic bounds to
% be set.  These default to the "Mercator" projection, centered
% at (longitude, latitude, rotation-angle) = (0, 0, 0) degrees,
% with no fixed bounds.
%    The "Transverse Mercator" projection can be computed
% from the "Mercator" mapping, using a center-of-projection
% of (LON, 0, 90), where LON is the desired central-longitude,
% in degrees.
%    The "Lambert Equal Area" projection is not conformal,
% but is useful for comparing the relative areas of adjacent
% grid cells.
 
% svn $Id$
%=======================================================================
% Copyright (C) 1999 Dr. Charles R. Denham, ZYDECO.
%  All Rights Reserved.
%   Disclosure without explicit written consent from the
%    copyright owner does not constitute publication.
%=======================================================================
 
% Version of 28-Dec-1999 11:00:24.
% Updated    28-Dec-1999 17:30:29.

helpdlg(help(mfilename), 'SeaGrid Projection')
