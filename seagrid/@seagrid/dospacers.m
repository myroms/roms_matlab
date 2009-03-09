function theResult = dospacers(self)

% seagrid/dospacers -- Compute the spacer effects.
%  dospacers(self) computes and draws the spacers of
%   the orthogonal grid for self, a "seagrid" object.
 
% svn $Id$
%=======================================================================
% Copyright (C) 1999 Dr. Charles R. Denham, ZYDECO.
%  All Rights Reserved.
%   Disclosure without explicit written consent from the
%    copyright owner does not constitute publication.
%=======================================================================
 
% Version of 15-Apr-1999 08:56:48.
% Updated    28-Jun-2000 17:06:00.

if nargout > 0, theResult = []; end
if nargin < 1, help(mfilename), return, end

theGridSize = psget(self, 'itsGridSize');
theSpacedGrids = psget(self, 'itsSpacedGrids');
theColor = psget(self, 'itsEdgeColor');
theEraseMode = psget(self, 'itsEraseMode');
theSpacerFlag = psget(self, 'itsSpacerFlag');
theSpacers = psget(self, 'itsSpacers');
theSpacerCounts = psget(self, 'itsSpacerCounts');
theSpacings = psget(self, 'itsSpacings');
theSpacedEdges = psget(self, 'itsSpacedEdges');
theSpacerMarker = psget(self, 'itsSpacerMarker');
theInterpFcn = psget(self, 'itsInterpFcn');
theInterpMethod = psget(self, 'itsInterpMethod');

theMarkerSize = get(0, 'DefaultLineMarkerSize');
theMarkerSize = 9;
theMarkerSize = 15;

switch theSpacerFlag
case 0
	theVisible = 'off';
otherwise
	theVisible = 'on';
end

if isempty(theSpacings)
	for i = 1:length(theGridSize)
		theSpacings{i} = linspace(0, 1, theGridSize(i)+1)
	end
	psset(self, 'itsSpacings', theSpacings)
end

TOLERANCE = sqrt(eps);

for k = 1:length(theGridSize)
	x = linspace(0, 1, length(theSpacings{k}));
	y = theSpacings{k};
	xi = linspace(0, 1, theSpacerCounts(k)+2);
	yi = feval(theInterpFcn, x, y, xi, theInterpMethod);
	for i = 2:length(yi)
		yi(i) = max(yi(i), yi(i-1)+TOLERANCE);
	end
	yi = yi / max(yi);
	dd{k} = yi;
end

psset(self, 'itsSpacerSpacings', dd)

dd1 = dd{1};
dd2 = dd{2};

f = findobj('Type', 'line', 'Tag', 'spacer');
if any(f), delete(f), end

TOLERANCE = sqrt(eps);

% The following gets confused when side #3
%  and/or side #4 is invoked.  I think that
%  we are mistakenly flip-flopping our spacers.

if ~isempty(theSpacedGrids)
	u = theSpacedGrids{1}; v = theSpacedGrids{2};
	u1 = u(:, 1); v1 = v(:, 1);
	u2 = u(end, :).'; v2 = v(end, :).';
	if theSpacedEdges(1) == 3
		u1 = u(:, end); v1 = v(:, end);
		u1 = flipud(u1); v1 = flipud(v1);
	end
	if theSpacedEdges(2) == 4
		u2 = u(1, :).'; v2 = v(1, :).';
		u2 = flipud(u2); v2 = flipud(v2);
	end
	z1 = u1 + v1*sqrt(-1);
	z2 = u2 + v2*sqrt(-1);
	d1 = cumsum([0; abs(diff(z1))]);
	for i = 2:length(d1)
		d1(i) = max(d1(i), d1(i-1)+TOLERANCE);
	end
	d1 = d1 / max(d1);
	d2 = cumsum([0; abs(diff(z2))]);
	for i = 2:length(d2)
		d2(i) = max(d2(i), d2(i-1)+TOLERANCE);
	end
	d2 = d2 / max(d2);
	zz1 = feval(theInterpFcn, d1, z1, dd1, theInterpMethod);
	zz2 = feval(theInterpFcn, d2, z2, dd2, theInterpMethod);
	uu1 = real(zz1); vv1 = imag(zz1);
	uu2 = real(zz2); vv2 = imag(zz2);
	theButtonDownFcn = 'psevent down';
	hold on
	for i = 1:length(uu1)
		h1(i) = plot(uu1(i), vv1(i), theSpacerMarker, ...
				'MarkerSize', theMarkerSize, ...
				'EraseMode', theEraseMode, 'Visible', theVisible, ...
				'ButtonDownFcn', theButtonDownFcn, 'Tag', 'spacer', ...
				'UserData', [uu1(i) vv1(i) dd1(i)]);
	end
	for i = 1:length(uu2)
		h2(i) = plot(uu2(i), vv2(i), theSpacerMarker, ...
				'MarkerSize', theMarkerSize, ...
				'EraseMode', theEraseMode, 'Visible', theVisible, ...
				'ButtonDownFcn', theButtonDownFcn, 'Tag', 'spacer', ...
				'UserData', [uu2(i) vv2(i) dd2(i)]);
	end
	set(h1([1 length(h1)]), 'Visible', 'off')
	set(h2([1 length(h2)]), 'Visible', 'off')
	hold off
	theSpacers = {h1, h2};
end

psset(self, 'itsSpacers', theSpacers)

if nargout > 0, theResult = self; end
