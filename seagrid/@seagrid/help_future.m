function help_future(self)

%    We intend for SeaGrid to be based solely on portable
% Matlab-5 m-files, independent of MathWorks toolboxes, and
% without any mex-files derived from lower-level languages.
% Presently, depending on the computer, three mex-files may
% be used, known as "mexrect", "mexsepeli", and "mexinside".
% The first can be replaced by "rect.m", and the third by
% "insidesafe.m".  The second can be inexactly mimicked by
% "fps.m", but an exact replacement will require additional
% work.
 
% svn $Id$
%=======================================================================
% Copyright (C) 1999 Dr. Charles R. Denham, ZYDECO.
%  All Rights Reserved.
%   Disclosure without explicit written consent from the
%    copyright owner does not constitute publication.
%=======================================================================
 
% Version of 25-May-1999 16:27:20.
% Updated    09-Jun-1999 08:10:11.

seagrid_helpdlg(help(mfilename), 'Future of SeaGrid Software')
