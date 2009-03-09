function help_spacing(self)

%    Small buttons called "spacers" along the edge of the SeaGrid
% boundary allow the grid-lines to be rearranged smoothly.  To see
% the spacers, use the "View/Spacers" menu item.  By dragging the
% spacers closer together, using mouse Button #1, the intervening
% grid lines themselves can be brought closer together.  To drag
% a spacer without forcing a full update of the window, use Button
% #2.  (Note: when lines are very closely spaced, they may disappear
% from the drawing, because of the "xor" drawing mode.  They do not,
% however, disappear from the computed grid.)
%    The number of spacers is controlled by the "View/Spacer Count"
% dialog.  Initially, there are five spacers along the two principal
% edges of the grid.
%    To restore the default spacing, use the "View/Default Spacing..."
% menu item to display the corresponding dialog.  The spacings
% along the two principal edges are entered as functions of "s",
% which uniformly spans the range of [0:1]. The dialog also permits
% the spacers to be moved to the opposite side of the grid.
 
% svn $Id$
%=======================================================================
% Copyright (C) 1999 Dr. Charles R. Denham, ZYDECO.
%  All Rights Reserved.
%   Disclosure without explicit written consent from the
%    copyright owner does not constitute publication.
%=======================================================================
 
% Version of 25-May-1999 16:27:20.
% Updated    04-Feb-2000 16:37:42.

seagrid_helpdlg(help(mfilename), 'SeaGrid Spacing')
