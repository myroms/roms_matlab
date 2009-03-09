function help_io(self)

%    Coastline data consist of counter-clockwise closed contours,
% separated by NaNs. The geographic positions must be stored in a
% Matlab mat-file, using the variable names "lon" and "lat", in degrees.
% The current coastline data can be replaced at any time, using the
% "SeaGrid/Load/Coastline" menu item.
%
%    Bathymetric data consist of arrays of points, stored in a mat-file
% under the names "xbathy", "ybathy", and "zbathy", corresponding to
% longitude, latitude (both degrees), and positive depth (arbitrary
% units), respectively.  The current bathymetric data can be replaced
% at any time, using "SeaGrid/Load/Bathymetry".
%
%    The current SeaGrid gridwork can be saved to a mat-file, using the
% "SeaGrid/Save" or "SeaGrid/Save As" menu item, then later retrieved
% with "Seagrid/Load" or "SeaGrid/Revert To Saved".  The grid and its
% associated data are stored in a Matlab "struct" named "s".  For an
% M-by-N grid, the grid-corners occupy an (M+1)-by-(N+1) array, while
% the gridded bathymetry and mask are stored as M-by-N arrays,
% corresponding to the grid-cell centers.
% 
%
%    To convert a SeaGrid file to Ecom or Scrum format, use the
% "seagrid2ecom" or "seagrid2scrum" routines.  See their respective
% "help" entries.
 
% svn $Id$
%=======================================================================
% Copyright (C) 1999 Dr. Charles R. Denham, ZYDECO.
%  All Rights Reserved.
%   Disclosure without explicit written consent from the
%    copyright owner does not constitute publication.
%=======================================================================
 
% Version of 25-May-1999 16:27:20.
% Updated    03-Apr-2000 11:50:02.

seagrid_helpdlg(help(mfilename), 'SeaGrid Input/Output')
