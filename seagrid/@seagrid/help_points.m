function help_points(self)

%    Use the "View" menu to highlight either the "Control Points"
% or the "Spacers".  The corner-points can be rolled to bring
% a particular starting corner (*) and its counter-clockwise
% neighbor into prominance.  The grid-line spacings can be
% reverted to uniform density along the two prominant
% edges, using the "View/Spacing Setup" menu item.
%
%    To add a draggable control-point, use button #1 on the
% grid boundary. To delete a control-point, use button #3.
%
%    Spacers along the grid boundary are used for fine-tuning
% the densities and positions of grid-lines.  Clicking button #1
% on a spacer increases the local grid-line density; clicking
% button #2 causes a decrease.  The adjustment factor can be
% set in the "View/Setup" dialog.  No grid lines are added.
% To move a grid-line manually within its neighborhood, drag
% with button #3.
 
% svn $Id$
%=======================================================================
% Copyright (C) 1999 Dr. Charles R. Denham, ZYDECO.
%  All Rights Reserved.
%   Disclosure without explicit written consent from the
%    copyright owner does not constitute publication.
%=======================================================================
 
% Version of 25-May-1999 16:27:20.
% Updated    27-Sep-1999 10:59:41.

seagrid_helpdlg(help(mfilename), 'SeaGrid Control Points')
