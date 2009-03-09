function help_start(self)

% [> PROJECTION <] -- Use "View/Setup" to select a projection.
% [> COASTLINE <] -- Use "SeaGrid/Load" to load a Mat-file containing
%      "lon" and "lat" variables. (Optional)
% [> BATHYMETRY <] -- Use "SeaGrid/Load" to load a Mat-file containing
%      "xbathy", "ybathy", and "zbathy" variables (longitude, latitude,
%      and positive depth, respectively).  (Optional)
% [> CORNERS <] -- Click with button #1 on four locations in counter
%      clockwise sequence to initialize the four corner-points of the
%      boundary.
% [> EDGES <] -- Click and drag with button #1 on any edge to insert
%      and move a new control-point.
% [> DRAG <] -- Drag with button #1 on a control-point to move it.
% [> DELETE <] -- Click button #3 on a control-point to delete it.
% [> DENSITY <] -- Use button #1 (button #2) on spacer-point to
%      increase (decrease) local grid cell density.  Drag with button
%      #3 to reposition a grid-line.
% [> ZOOM <] -- Click elsewhere to invoke "zoomsafe" action.
 
% seagrid/help_start

% svn $Id$
%=======================================================================
% Copyright (C) 1999 Dr. Charles R. Denham, ZYDECO.
%  All Rights Reserved.
%   Disclosure without explicit written consent from the
%    copyright owner does not constitute publication.
%=======================================================================
 
% Version of 16-Apr-1999 07:40:25.
% Updated    28-Dec-1999 18:11:56.

seagrid_helpdlg(help(mfilename), 'Getting Started With SeaGrid')
