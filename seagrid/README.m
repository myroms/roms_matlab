function README

help(mfilename)

%
% ************* Seagrid Tutorial *************
%
%    Start by executing "seagrid" at the Matlab prompt.
% Its many features are described in the "help" entries,
% accessible through the "Seagrid/Help" menu.
%
%    Select a projection by invoking the "View/Setup"
% menu item to display the "SeaGrid Setup" dialog.
% The default is "Mercator".
%
%    The next task is to open a coastline file, if any,
% followed by a bathymetry file, if any.  The required
% formats are given in the "help".  These data will be
% plotted on the current projection.  The plot units
% can be toggled between "Degrees" and "Kilometers"
% by using the "View/Units" menu item.
%
%    Then, four corner-points need to be placed on the map
% by simple clicking.  These points are draggable, so their
% exact placement is not critical.
%
%    Macintosh clicking: #1 = simple click; #2 = shift-click;
% and #3 = option-click.
%
%    Once the four corner-points have been placed, the map
% becomes zoomable, using the "zoomsafe" protocol.  See
% "help zoomsafe".
%
%    To add an edge-point, click on any edge, dragging if
% desired.  As with corner-points, edge-points are draggable
% at any time.  Also, an edge-point can be deleted by clicking
% on it with button #3.  The boundary is displayed with cubic
% splines that pass through the corner-points and edge-points.
%
%    The number of gridlines is adjustable through the "Setup"
% dialog (select the "View/Setup" menu item).
%
%    The grid cell spacing can be manipulated by invoking
% the "View/Spacers" menu item.  Clicking (button #1) and
% dragging a spacer will rearrange the grid lines smoothly
% to mirror how closely spaced the spacers themselves are.
% The number of spacers can be changed with the
% "View/Spacer Count" menu item.
%
%    The default grid cell spacing can be expressed as a
% function through the "View/Default Spacing" menu item.
% The subsequent dialog accepts functions of "s", the
% relative distance [0:1] along each of the two principal
% edges of the boundary.  The beginning of the first
% principal edge is marked by "*", and the boundary always
% runs counter-clockwise always.  This geometry can be rolled
% with the "View/Roll" menu item.
%
%    The current grid can be saved for later retrieval, using
% the "Seagrid/Save" or "Seagrid/Save As" menu items.  The
% most recently saved grid is available through the
% "SeaGrid/Revert To Saved" menu item.
%
%    The grid can be masked with the "Compute/Land Mask" menu,
% to mark those cells whose centers fall on land.  Also, the
% "Mask Tool" can be invoked through the "Compute/Mask Tool"
% menu item to allow manual toggling of the mask.  Finally,
% the "Compute/Erase Mask" menu can be used to erase the
% current mask.
%
%   NOTE: the manual masking changes do not become effective until
% the "Mask Tool" is dismissed by choosing the "Compute/Mask Tool"
% menu item again.
%
%   To convert the saved SeaGrid output to a SCRUM or ECOM file,
% execute "seagrid2scrum" or "seagrid2ecom".  Both routines will
% prompt for filenames if none are provided on the command-line.
% See "help" for each.

% svn $Id$
%=======================================================================
% Copyright (C) 1999 Dr. Charles R. Denham, ZYDECO.
%  All Rights Reserved.
%   Disclosure without explicit written consent from the
%    copyright owner does not constitute publication.
%=======================================================================
%
% Version of 05-Aug-1999 14:26:12.
% Updated    23-Jul-2001 13:46:26.

help(mfilename)
