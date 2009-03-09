function theResult = domask(self, needsUpdate)

% seagrid/domask -- Compute the "seagrid" mask.
%  domask(self) computes the mask for self,
%   a "seagrid" object.  Non-zero mask-values
%   denote land; zeros denotes water.
 
% svn $Id$
%=======================================================================
% Copyright (C) 1999 Dr. Charles R. Denham, ZYDECO.
%  All Rights Reserved.
%   Disclosure without explicit written consent from the
%    copyright owner does not constitute publication.
%=======================================================================
 
% Version of 10-May-1999 09:40:20.
% Updated    14-Nov-2000 11:24:14.

if nargout > 0, theResult = self; end
if nargin < 1, help(mfilename), return, end
if nargin < 2, needsUpdate = psget(self, 'itsNeedsUpdate'); end

theGridSize = psget(self, 'itsGridSize');

theCoastline = findobj(gcf, 'Tag', 'coastline');
if isempty(theCoastline), return, end
cx = get(theCoastline, 'XData');
cy = get(theCoastline, 'YData');

theSpacedGrids = psget(self, 'itsSpacedGrids');
u = theSpacedGrids{1};
v = theSpacedGrids{2};
u = interp2(u, 1);
v = interp2(v, 1);
utemp = u(2:2:end-1, 2:2:end-1);  % Grid-centers.
vtemp = v(2:2:end-1, 2:2:end-1);

useMexFile = psget(self, 'itUsesMexFile');

hasMex = (exist('mexinside', 'file') == 3);
theFcn = 'insidesafe';
if hasMex & useMexFile, theFcn = 'mexinside'; end

if needsUpdate
	busy
	f = find(~isfinite(cx) | ~isfinite(cy));
	f = f(:).';
	if ~any(f), f = [0 length(cx)+1]; end
	if f(1) ~= 1, f = [0 f]; end
	if f(end) ~= length(cx), f(end+1) = length(cx)+1; end
	if isunix
		shouldAlert = (prod(theGridSize) > 1000);
	else
		shouldAlert = (prod(theGridSize) > 100);
	end
	if shouldAlert
		h = warndlg('Please wait ...', 'Computing Mask');
		drawnow
	end
	theMask = zeros(size(utemp));
	for i = 2:length(f)
		g = find(theMask == 0);
		if ~any(g), break, end
		j = f(i-1)+1:f(i)-1;
		if length(j) > 2
			theMask(g) = feval(theFcn, utemp(g), vtemp(g), cx(j), cy(j));
		end
	end
	if shouldAlert & ishandle(h), delete(h), end
	
	theLand = ~~theMask;
	theWater = ~theLand;
	
	psset(self, 'itsMask', theMask)
	psset(self, 'itsWater', theWater)
	psset(self, 'itsLand', theLand)
	idle
else
	theMask = psget(self, 'itsMask');
end

theFigure = ps(self);
setsafe(0, 'CurrentFigure', theFigure)

if nargout > 0, theResult = self; end
