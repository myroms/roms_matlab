function theResult = dospacings(self)

% seagrid/dospacings -- Compute spacings from spacer locations.
%  dospacings(self) computes the grid spacings from the spacer
%   locations, on behalf of self, a "seagrid" object.  The
%   density of grid-lines in the vicinity of the selected
%   point is doubled (click) or halved (shift-click).
 
% svn $Id$
%=======================================================================
% Copyright (C) 1999 Dr. Charles R. Denham, ZYDECO.
%  All Rights Reserved.
%   Disclosure without explicit written consent from the
%    copyright owner does not constitute publication.
%=======================================================================
 
% Version of 23-Apr-1999 12:51:23.
% Updated    11-Feb-2000 16:45:12.

if nargout > 0, theResult = self; end
if nargin < 1, help(mfilename), return, end

theSpacers = psget(self, 'itsSpacers');
if isempty(theSpacers), return, end

sel = psget(self, 'SelectionType');

theGridSize = psget(self, 'itsGridSize');
theInterpFcn = psget(self, 'itsInterpFcn');
theInterpMethod = psget(self, 'itsInterpMethod');
theSpacerSpacings = psget(self, 'itsSpacerSpacings');
theSpacerCounts = psget(self, 'itsSpacerCounts');
theSpacings = psget(self, 'itsSpacings');
theSpacedEdges = psget(self, 'itsSpacedEdges');

if isempty(theSpacings)
	theSpacings = cell(1, length(theGridSize));
	for k = 1:length(theSpacings)
		theSpacings{k} = linspace(0, 1, theGridSize(k)+1);
	end
	psset(self, 'itsSpacings', theSpacings)
end

if isempty(theSpacedEdges)
	theSpacedEdges = [1 2];
	psset(self, 'itsSpacedEdges', theSpacedEdges)
end

for k = 1:length(theGridSize)
	h = theSpacers{k};
	for i = 1:length(h)
		if ~ishandle(h(i))
			psset(self, 'itsSpacers', [])
			psset(self, 'itsSpacings', [])
			return
		end
	end
end

for k = 1:length(theGridSize)
	h = theSpacers{k};
	if isempty(h)
		psset(self, 'itsSpacings', [])
	end
	selected = 0;
	for i = 1:length(h)
		if isequal(get(h(i), 'Selected'), 'on')
			selected = i;
			break;
		end
	end
	if selected > 1 & selected < length(h)
		s = theSpacerSpacings{k};
		s = s(:).';   % Make row-vector.
		switch sel
		case {'normal', 'extend'}   % Drag a spacer.
			z = zeros(size(h));
			zold = zeros(size(h));
			didMove = 0;
			for i = 1:length(h)
				x = get(h(i), 'XData');   % New x location.
				y = get(h(i), 'YData');   % New y location.
				udata = get(h(i), 'UserData');   % Old [x y] location.
				xyold = udata(1:2);
				if any(xyold(1:2) ~= [x y]) & i == selected
					didMove = 1;
				end
				z(i) = x + sqrt(-1)*y;
				zold(i) = xyold(1) + sqrt(-1)*xyold(2);
				setsafe(h(i), 'UserData', [x y])
			end
			if didMove
				d = abs(diff(z));
				a = d(selected-1);
				b = d(selected);
				frac = a ./ (a + b);
				s(selected) = (1-frac) .* s(selected-1) + frac * s(selected+1);
				udata = get(h(selected), 'UserData');
				udata = [udata s(selected)];
				setsafe(h(selected), 'UserData', udata)
			end
		end
		theSpacerSpacings{k} = s;
		x = linspace(0, 1, length(s));
		y = s;
		xi = linspace(0, 1, theGridSize(k)+1);
		yi = splinesafe(x, y, xi, 1);
		for i = 2:length(yi)
			yi(i) = max(yi(i), yi(i-1));
		end
		theSpacings{k} = yi;
	end
end

% Interpolate if grid-size has changed.

TOLERANCE = sqrt(eps);

for k = 1:length(theSpacings)
	if length(theSpacings{k}) ~= theGridSize(k)+1
		x = linspace(0, 1, length(theSpacings{k}));
		y = theSpacings{k};
		xi = linspace(0, 1, theGridSize(k)+1);
		yi = feval(theInterpFcn, x, y, xi, theInterpMethod);
		for i = 2:length(yi)
			yi(i) = max(yi(i), yi(i-1)+TOLERANCE);
		end
		yi = yi / max(yi);
		theSpacings{k} = yi;
	end
end

psset(self, 'itsSpacings', theSpacings)
psset(self, 'itsSpacerSpacings', theSpacerSpacings)

if nargout > 0, theResult = self; end
