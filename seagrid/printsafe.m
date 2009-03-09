function PrintSafe(theFigures, theOptions)

% PrintSafe -- Safely print figures.
%  PrintSafe(theFigure, 'theOptions') prints theFigures
%   (default = current-figure), using 'theOptions'.  This
%   routine temporarily disables the 'ResizeFcn' property
%   and applies the '-noui' option automatically.  Any
%   additional options are provided in a string that is
%   appended to the "print" command.  For example, use
%   "printsafe(gcf, '-loose') to impose the 'PaperPosition'
%   property of each figure onto the corresponding PostScript
%   "BoundingBox".
%  PrintSafe('theOptions') performs "printsafe(gcf, theOptions)".

% svn $Id$
%=======================================================================
% Copyright (C) 1999 Dr. Charles R. Denham, ZYDECO.
%  All Rights Reserved.
%   Disclosure without explicit written consent from the
%    copyright owner does not constitute publication.
%=======================================================================

if nargin < 1, help(mfilename), theFigures = gcf; end
if nargin < 2, theOptions = ''; end

if isstr(theFigures)
   if theFigures(1) ~= '-'
      theFigures = eval(theFigures);
     else
      theOptions = theFigures;
      theFigures = gcf;
   end
end

f = findobj('Type', 'figure');
if ~any(f), return, end

theGCF = gcf;

for k = 1:length(theFigures)
   figure(theFigures(k))
   theResizeFcn = get(theFigures(k), 'ResizeFcn');
   set(theFigures(k), 'Resizefcn', '')
   theCommand = ['print -noui -f' num2str(theFigures(k)) ' ' theOptions];
   eval(theCommand, ['disp(['' ## Not successful '' theCommand])'])
   set(theFigures(k), 'Resizefcn', theResizeFcn)
end

figure(theGCF)
