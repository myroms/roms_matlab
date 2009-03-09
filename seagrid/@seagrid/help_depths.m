function help_depths(self)

%    Select the "Compute/Grid Depths" menu item to compute and display
% depths for the grid-cell centers from the current bathymetry
% and coastline data.  Cells on land, plus those for which no
% depth can be determined, are displayed transparently.
% The state of the mask can be toggled with mouse button #1.
% While the "Compute/Grid Depths" item is active, the depths are
% recomputed automatically whenever the grid geometry changes.
 
% svn $Id$
%=======================================================================
% Copyright (C) 1999 Dr. Charles R. Denham, ZYDECO.
%  All Rights Reserved.
%   Disclosure without explicit written consent from the
%    copyright owner does not constitute publication.
%=======================================================================
 
% Version of 17-Jun-1999 16:42:56.
% Updated    06-Jul-2000 14:17:19.

seagrid_helpdlg(help(mfilename), 'SeaGrid Depths')
