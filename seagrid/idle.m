function Idle(theFigure)

% Idle -- Set the arrow-cursor.
%  Idle(theFigure) sets the arrow-cursor in theFigure.
%   The companion routine is "Busy".
 
% svn $Id$
%=======================================================================
% Copyright (C) 1996 Dr. Charles R. Denham, ZYDECO.
%  All Rights Reserved.
%   Disclosure without explicit written consent from the
%    copyright owner does not constitute publication.
%=======================================================================

if ~any(findobj('Type', 'figure')), return, end

if nargin < 1, theFigure = gcf; end

set(theFigure, 'Pointer', 'arrow');
