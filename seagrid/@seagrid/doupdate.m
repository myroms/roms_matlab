function theResult = doupdate(self, needsUpdate)

% seagrid/doupdate -- Update a "seagrid" object.
%  doupdate(self) updates self, a "seagrid" object.
%   If the "itNeedsUpdate" flag of self is logically
%   true, the grid is recomputed as part of the
%   updating protocol.
%  update(self, needsUpdate) uses the given flag
%   to override the 'itNeedsUpdate' flag in self.
 
% svn $Id$
%=======================================================================
% Copyright (C) 1999 Dr. Charles R. Denham, ZYDECO.
%  All Rights Reserved.
%   Disclosure without explicit written consent from the
%    copyright owner does not constitute publication.
%=======================================================================
 
% Version of 07-Apr-1999 21:19:19.
% Updated    23-Apr-2001 16:00:02.

theProjection = psget(self, 'itsProjection');
theProjectionCenter = psget(self, 'itsProjectionCenter');
sg_proj(theProjection, theProjectionCenter)

theFigure = ps(self);
setsafe(0, 'CurrentFigure', theFigure)

theAxes = findobj(theFigure, 'Type', 'axes', 'Tag', 'seagrid');
if isempty(theAxes)
	theAxes = gca;
end
setsafe(theAxes, 'Tag', 'seagrid')
setsafe(theFigure, 'CurrentAxes', theAxes)

theTitle = get(theAxes, 'Title');
setsafe(theTitle, 'String', [theProjection ' Projection'])

thePoints = psget(self, 'itsPoints');
if length(thePoints) < 4
	doticks(self)
	setsafe(theFigure, 'Pointer', 'crosshair')
	if nargout > 0, theResult = self; end
	return
end

busy

if nargin < 2
	needsUpdate = psget(self, 'itNeedsUpdate');
else
	psset(self, 'itNeedsUpdate', needsUpdate)
end

if verbose(self)
	disp([' ## ' mfilename ' ' int2str(needsUpdate)])
end

theOrthogonalityFlag = psget(self, 'itsOrthogonalityFlag');
theBathymetryFlag = psget(self, 'itsBathymetryFlag');
theMaskingFlag = psget(self, 'itsMaskingFlag');
theMaskTool = psget(self, 'itsMaskTool');
theMaskToolFlag = psget(self, 'itsMaskToolFlag');

if any(theMaskTool)
	domasktool(self)   % Turn it off.
end

h = findobj(theFigure, 'Tag', 'masktool');
if any(h), delete(h), end

h = findobj(theFigure, 'Tag', 'grid2mask');
if any(h), delete(h), end

psset(self, 'itsMaskTool', []);
theMaskTool = psget(self, 'itsMaskTool');

dospacings(self);   % <== Step #1 <==

h = findobj(theFigure, 'Tag', 'orthogonality');
if any(h), delete(h), end
h = findobj(theFigure, 'Tag', 'gridded-bathymetry');
if any(h), delete(h), end
h = findobj(theFigure, 'Tag', 'contoured-bathymetry');
if any(h), delete(h), end
h = findobj(theFigure, 'Tag', 'mask');
if any(h), delete(h), end
h = findobj(theFigure, 'Type', 'line', 'Tag', 'grid-line');
if any(h), delete(h), end
h = findobj(theFigure, 'Type', 'line', 'Tag', 'edge');
if any(h), delete(h), end
h = findobj(theFigure, 'Type', 'line', 'Tag', 'spacer');
if any(h), delete(h), end

dogrid(self, needsUpdate);   % <== Step #2 <==

f = findobj(theFigure, 'Tag', 'Colorbar');
if any(f), delete(f), end

theAxes = gca;

theColorBarFlag = 0;
if any(theOrthogonalityFlag)
	self = doorthogonality(self);   % <== Step #3 <==
	theColorBarFlag = 0;   % Does its own colorbar.
else
	if any(theMaskingFlag)
		needsMask = 1;
		self = domask(self, needsMask);   % <== Step #4 <==
	end
	if any(theBathymetryFlag)
		needsBathymetry = 1;
		self = dobathymetry(self, needsBathymetry);   % <== Step #5 <==
	end
	if any(theMaskToolFlag)
		theMaskTool = psget(self, 'itsMaskTool');
		if any(theMaskTool)
			domasktool(self)   % Delete it.
		end
		self = domasktool(self);   % <== Step #6 <==
		theColorBarFlag = 1;
		psset(self, 'itsMaskingFlag', 0)
		psset(self, 'itsBathymetryFlag', 0)
	end
end

setsafe(0, 'CurrentFigure', theFigure)

doedges(self);   % <== Step #7 <==
dospacers(self);   % <== Step #8 <==
dopoints(self);   % <== Step #9 <==
doticks(self);   % <== Step #10 <==

% We have to activate the colorbar here
%  whenever we use the masktool.

if any(theColorBarFlag)
	setsafe(gca, 'CLimMode', 'auto')
	colorbar
end

if any(theColorBarFlag) & any(theMaskToolFlag) & ...
		~any(theOrthogonalityFlag)
	h = colorbar;
	setsafe(get(h, 'Ylabel'), 'String', 'Depth')
end

axes(theAxes)

theWBDF = ...
	'if zoomsafe(''down''), doticks(ps(gcf), 1), end';
	
setsafe(theFigure, 'WindowButtonDownFcn', theWBDF)

needsUpdate = 0;
psset(self, 'itNeedsUpdate', needsUpdate)

idle

if nargout > 0, theResult = self; end
