function help_about(self, varargin)

%    SeaGrid provides a graphical user interface for constructing
% a one-layer, orthogonal, curvilinear, geographic grid, which
% is flat on the top and bathymetrically controlled on the bottom.
%
%    Information on various SeaGrid topics is available in the
% "Help" menu.
%
%    Designed and programmed by Chuck Denham, with the assistance
% of Rich Signell (author of "cgridgen"), John Evans (author of
% "gridgen"), Courtney Harris, Rob Hetland, Jason Hyatt.

% svn $Id$
%=======================================================================
% Copyright (C) 1999 Dr. Charles R. Denham, ZYDECO.
%  All Rights Reserved.
%   Disclosure without explicit written consent from the
%    copyright owner does not constitute publication.
%=======================================================================

% Version of 03-Jun-1999 22:20:01.
% Updated    23-Dec-1999 09:19:25.

seagrid_helpdlg(help(mfilename), 'About SeaGrid')
