function help_setup_dialog(self)

%    The SeaGrid "Setup" dialog can be displayed by
% selecting the "View/Setup..." menu item.  The dialog
% offers selectors for the map projection (default =
% "Mercator"); geographic bounds (optional; 
% defaults = []); the number of cells along edges #1
% and #2 of the grid (defaults = 10); the depth range
% (defaults = [0 Inf]); the end-condition for the cubic
% splines that are used throughout SeaGrid (default =
% "end-slope flag on"); and whether or not to invoke
% Mex-files (default = "yes").
%    SeaGrid always performs a full update whenever
% the dialog is accepted via its "Okay" button.  No
% update is performed if the dialog is dismissed via
% its 'Cancel" button or "go-away" box.
 
% svn $Id$
%=======================================================================
% Copyright (C) 1999 Dr. Charles R. Denham, ZYDECO.
%  All Rights Reserved.
%   Disclosure without explicit written consent from the
%    copyright owner does not constitute publication.
%=======================================================================
 
% Version of 27-Dec-1999 06:42:00.
% Updated    27-Dec-1999 06:42:00.

helpdlg(help(mfilename), 'SeaGrid Setup Dialog')
