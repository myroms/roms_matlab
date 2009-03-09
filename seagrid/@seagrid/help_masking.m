function help_masking(self)

%    Select the "Compute/Land Mask" menu item to compute
% and display the land-mask for the grid-cell centers, using
% the current coastline data.  Cells on land are displayed
% transparently.  The state of the mask can be toggled with
% mouse button #1.  While the "Compute/Land Mask" item is
% active, the mask is recomputed automatically whenever
% the grid geometry changes.
 
% svn $Id$
%=======================================================================
% Copyright (C) 1999 Dr. Charles R. Denham, ZYDECO.
%  All Rights Reserved.
%   Disclosure without explicit written consent from the
%    copyright owner does not constitute publication.
%=======================================================================
 
% Version of 17-Jun-1999 16:42:56.
% Updated    17-Aug-2000 10:53:59.

seagrid_helpdlg(help(mfilename), 'SeaGrid Masking')
