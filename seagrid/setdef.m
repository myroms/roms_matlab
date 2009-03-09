function setdef(theFunction)

% setdef -- Switch to directory of a function.
%  setdef(theFunction) switches to the directory
%   that contains theFunction.
%  setdef (no argument) displays pwd.

% svn $Id$
%=======================================================================
% Copyright (C) 1996 Dr. Charles R. Denham, ZYDECO.
%  All Rights Reserved.
%   Disclosure without written consent from the
%    copyright owner does not constitute publication.
%=======================================================================

if nargin > 0
   w = which(theFunction);
   f = findstr(w, filesep);
   if any(f)
      w = w(1:f(length(f)));
      eval(['cd(''' w(1:f(length(f))) ''')'])
   end
end
disp([' ## ' pwd])
