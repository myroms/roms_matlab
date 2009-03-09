function help_math(self)

%    Starting with an irregular physical boundary on a conformal
% (angle-preserving) projection, such as Mercator, SeaGrid maps it
% to a rectangle, using the conformal "wire-straightening" scheme
% of Ives and Zacharias (1987).  The X coordinates of the former
% are then distributed around the perimeter of the latter.  The
% interior of the rectangle is filled by a fast Poisson solver.
% The process is repeated independently for the Y coordinate.
%    The bathymetry is gridded at the cell-centers, using the Matlab
% "griddata" routine.  The mask is computed at the cell-centers: a
%  cell is considered to be fully on-land if its center is on-land.
%
%    Reference: David C. Ives and Robert M. Zacharias, Conformal
% Mapping and Orthogonal Grid Generation (Paper No. AA-87-2057),
% AIAA/SAE/ASME/ASEE 23rd Joint Propulsion Conference, San Diego,
% California, June 29-July 2,1987.
 
% svn $Id$
%=======================================================================
% Copyright (C) 1999 Dr. Charles R. Denham, ZYDECO.
%  All Rights Reserved.
%   Disclosure without explicit written consent from the
%    copyright owner does not constitute publication.
%=======================================================================
 
% Version of 25-May-1999 16:27:20.
% Updated    09-Jun-1999 08:53:39.

seagrid_helpdlg(help(mfilename), 'SeaGrid Mathematics')
