function theResult = dogrid(self, needsUpdate)

% seagrid/dogrid -- Compute the orthogonal grid.
%  dogrid(self, needsUpdate) computes and draws
%   the curvilinear orthogonal grid, using the
%   control-points of self, a "seagrid" object.
 
% svn $Id$
%=======================================================================
% Copyright (C) 1999 Dr. Charles R. Denham, ZYDECO.
%  All Rights Reserved.
%   Disclosure without explicit written consent from the
%    copyright owner does not constitute publication.
%=======================================================================
 
% Version of 15-Apr-1999 08:56:48.
% Updated    21-Sep-2000 15:03:33.

if nargout > 0, theResult = self; end

if nargin < 1, help(mfilename), return, end
if nargin < 2, needsUpdate = psget(self, 'itNeedsUpdate'); end

if nargout > 0, theResult = self; end

theEraseMode = psget(self, 'itsEraseMode');
theGridSize = psget(self, 'itsGridSize');
theSpacings = psget(self, 'itsSpacings');

u = []; v = [];
uu = []; vv = [];

theGrids = psget(self, 'itsGrids');
if ~isempty(theGrids)
	u = theGrids{1}; v = theGrids{2};
end

theSpacedGrids = psget(self, 'itsSpacedGrids');
if ~isempty(theSpacedGrids)
	uu = theSpacedGrids{1}; vv = theSpacedGrids{2};
end

theInterpFcn = psget(self, 'itsInterpFcn');
theInterpMethod = psget(self, 'itsInterpMethod');

c = computer;
isMac = isequal(lower(c(1:3)), 'mac');
doUseMexFile = psget(self, 'itUsesMexFile');
doUseMexFile = (doUseMexFile & ~isMac);

if isempty(u), needsUpdate = 1; end
	
if needsUpdate
	
	busy
	
	if isunix
		shouldAlert = (prod(theGridSize) > 10000);
	else
		shouldAlert = (prod(theGridSize) > 1000);
	end
	if shouldAlert
		h = warndlg('Please wait ...', 'Computing Grid');
		drawnow
	end
	
	gridcalc(self)   % <== Calculate grid.
	
	if shouldAlert & ishandle(h), delete(h), end
	
	idle
end

theGrids = psget(self, 'itsGrids');

if ~isempty(theGrids)
	u = theGrids{1}; v = theGrids{2};
end

theSpacedEdges = psget(self, 'itsSpacedEdges');
if isempty(theSpacedEdges)
	theSpacedEdges = [1 2];
	psset(self, 'itsSpacedEdges', theSpacedEdges)
end

if isempty(u)
	psset(self, 'itsGrids', [])
	h = warndlg('Please rearrange control-points.', 'Bad Shape');
	pause(2)
	if ishandle(h), delete(h), end
	if nargout > 0, theResult = self; end
	return
end

% Apply spacings for the plotted grid.
%  The spacings are the relative tic-locations
%  [0:1] on the computed grid where we want to
%  interpolate and show grid lines.  The spacings
%  are defined by the positions of the spacers
%  in the "seagrid" window.

uu = u;
vv = v;
[m, n] = size(u);

% Interpolate to even-spacing along the
%  curved paths of the two reference sides.
%  N.B. -- Needs to be changed to "spline"
%  interpolation.

if ~doUseMexFile

	sm = linspace(0, 1, m);
	sn = linspace(0, 1, n);

	if ~isempty(u) & 1
		zz = uu + sqrt(-1)*vv;
		switch theSpacedEdges(1)
		case 1
			s1 = rescale([0; cumsum(abs(diff(zz(:, 1))))]);
		case 3
			s1 = rescale([0; cumsum(abs(diff(zz(:, end))))]);
			s1 = 1 - s1(end:-1:1);
		end
		switch theSpacedEdges(2)
		case 2
			s2 = rescale([0 cumsum(abs(diff(zz(end, :))))]).';
		case 4
			s2 = rescale([0 cumsum(abs(diff(zz(1, :))))]).';
			s2 = 1 - s2(end:-1:1);
		end
		uu = feval(theInterpFcn, s1, uu, sm, theInterpMethod).';   % Down.
		vv = feval(theInterpFcn, s1, vv, sm, theInterpMethod).';
		uu = feval(theInterpFcn, s2, uu, sn, theInterpMethod).';   % Across.
		vv = feval(theInterpFcn, s2, vv, sn, theInterpMethod).';
	end

% Interpolate to manually-adjusted spacings.

	if ~isempty(u) & 1
		if isempty(theSpacings)
			theSpacings = {sm, sn};   % <== Note.
		end
		s1 = rescale(theSpacings{1});
		s2 = rescale(theSpacings{2});
		if theSpacedEdges(1) == 3
			s1 = 1 - s1(end:-1:1);
		end
		if theSpacedEdges(2) == 4
			s2 = 1 - s2(end:-1:1);
		end
		ss1 = linspace(0, 1, length(s1));
		ss2 = linspace(0, 1, length(s2));
		s1 = feval(theInterpFcn, ss1, s1, sm, theInterpMethod);
		s2 = feval(theInterpFcn, ss2, s2, sn, theInterpMethod);
		theSpacings = {s1, s2};   % <== Note.
		uu = feval(theInterpFcn, sm, uu, s1, theInterpMethod).';   % Down.
		vv = feval(theInterpFcn, sm, vv, s1, theInterpMethod).';
		uu = feval(theInterpFcn, sn, uu, s2, theInterpMethod).';   % Across.
		vv = feval(theInterpFcn, sn, vv, s2, theInterpMethod).';
	end

end

psset(self, 'itsSpacedGrids', {uu vv})
% psset(self, 'itsSpacings', {theSpacings})   % <== Note.

% Plot the grid.

if ~isempty(uu)
	[m, n] = size(uu);
	u1 = uu(:, 2:end-1); v1 = vv(:, 2:end-1);
	u2 = uu(2:end-1, :).'; v2 = vv(2:end-1, :).';
	hold on
	oldGridLines = findobj('Type', 'line', 'Tag', 'grid-line');
	if any(oldGridLines), delete(oldGridLines), end
	theGridLines{1} = plot(u1, v1, 'g-', ...
			'EraseMode', theEraseMode, 'Tag', 'grid-line');
	theGridLines{2} = plot(u2, v2, 'g-', ...
			'EraseMode', theEraseMode, 'Tag', 'grid-line');
	hold off
end

if nargout > 0, theResult = self; end
