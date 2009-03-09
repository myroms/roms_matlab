function help_bugs(self)

% 1. Load: Coastline input must preceed bathymetry.
% 2. Color: Colors may clash, due to the "xor" drawing
%       mode used in SeaGrid.  Selecting "View/Refresh"
%       may improve the appearance.
% 3. Mask: An empty mask may not show its grid-lines.
% 4. Pixels: Dots occasionally get orphaned.  Use
%       "View/Refresh" to redraw.
% 5. Macintosh: Grids are interpolated from a Fourier-
%       Poisson solution for uniform grid-lines, causing
%       some grid crossings to deviate from othogonality.
% 6. Unix: Matlab GUI behaviors are buggy in Unix.  The
%       "rubber-rectangle" will be invisible when doing
%       zooming with a click/drag motion.  The "tab" key
%       will not work for walking through dialogs.  Menu
%       accelerator-keys often fail to perform.
% 7. WinDoze: Modal dialogs sometimes do not lay-out
%       cleanly.

% svn $Id$
%=======================================================================
% Copyright (C) 1999 Dr. Charles R. Denham, ZYDECO.
%  All Rights Reserved.
%   Disclosure without explicit written consent from the
%    copyright owner does not constitute publication.
%=======================================================================

% Version of 03-Jun-1999 22:20:01.
% Updated    07-Mar-2000 22:54:26.

seagrid_helpdlg(help(mfilename), 'Known Bugs/Features In SeaGrid')
