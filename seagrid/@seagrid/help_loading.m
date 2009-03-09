function help_loading(self)

%    Coastline data consist of counter-clockwise closed contours,
% separated by NaNs. The geographic positions must be stored in a
% Matlab mat-file, using the variable names "lon" and "lat", in degrees.
% The current coastline data can be replaced at any time, using the
% "SeaGrid/Load/Coastline" menu item.
%
%    Bathymetric data consist of arrays of points, stored as three
% columns in an ASCII file, or stored in a mat-file as variables
% "xbathy", "ybathy", and "zbathy", corresponding to longitude,
% latitude (both degrees), and positive depth (arbitrary units),
% respectively.  The current bathymetric data can be replaced
% at any time, using "SeaGrid/Load/Bathymetry".
%
%    Boundary data can be loaded from a three-column ASCII file of
% of [lon lat isCorner], where isCorner is 1 for each corner and
% 0 for all others.  A corner-point must be the first element in
% the file.  Do NOT repeat the first point at the end.

% svn $Id$
%=======================================================================
% Copyright (C) 1999 Dr. Charles R. Denham, ZYDECO.
%  All Rights Reserved.
%   Disclosure without explicit written consent from the
%    copyright owner does not constitute publication.
%=======================================================================
 
% Version of 25-May-1999 16:27:20.
% Updated    04-Apr-2000 11:10:41.

seagrid_helpdlg(help(mfilename), 'SeaGrid Loading')
