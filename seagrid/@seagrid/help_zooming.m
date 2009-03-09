function help_zooming(self)

%    The SeaGrid window uses "zoomsafe" for zooming and
% panning.  Click in white-space or on inactive objects
% to zoom.  To invoke a rubber-rectangle, drag while the
% mouse button is down.  Existing callbacks in graphical
% objects take precedence.  Zooming temporarily overrides 
% any geographic bounds that have been set manually.
%
%   Button #1      Zoom in, centered on click or rectangle.
%   Button #2      Zoom out, centered on click or rectangle.
%   Button #3      Center on click or rectangle without zooming.
%   Double-click   Zoom fully out.
%
%   See "help zoomsafe" for more information.
 
% svn $Id$
%=======================================================================
% Copyright (C) 1999 Dr. Charles R. Denham, ZYDECO.
%  All Rights Reserved.
%   Disclosure without explicit written consent from the
%    copyright owner does not constitute publication.
%=======================================================================
 
% Version of 25-May-1999 16:27:20.
% Updated    27-Dec-1999 07:04:38.

seagrid_helpdlg(help(mfilename), 'SeaGrid Zooming')
