function theResult = gridcalc(self)

% gridcalc -- Calculate the basic grid for "seagrid".
%  gridcalc(self) calculates the grid for self, a
%   "seagrid" object.
 
% svn $Id$
%=======================================================================
% Copyright (C) 1999 Dr. Charles R. Denham, ZYDECO.
%  All Rights Reserved.
%   Disclosure without explicit written consent from the
%    copyright owner does not constitute publication.
%=======================================================================
 
% Version of 19-May-1999 05:55:33.
% Updated    14-Oct-1999 16:38:11.

if nargout > 0, theResult = self; end
if nargin < 1, help(mfilename), return, end

thePoints = psget(self, 'itsPoints');
theCornerTag = psget(self, 'itsCornerTag');
theGridSize = psget(self, 'itsGridSize');
theEndSlopeFlag = psget(self, 'itsEndSlopeFlag');
theInterpFcn = psget(self, 'itsInterpFcn');
theInterpMethod = psget(self, 'itsInterpMethod');
useMexFile = psget(self, 'itUsesMexFile');

% Check for Macintosh; no mex-files yet.

c = computer;
isMac = isequal(lower(c(1:3)), 'mac');
useMexFile = (useMexFile & ~isMac);

x = zeros(1, length(thePoints));
y = zeros(1, length(thePoints));

if min(size(thePoints)) == 1
	theCorners = [];
	for k = 1:length(thePoints)
		h = thePoints(k);
		x(k) = get(h, 'XData');
		y(k) = get(h, 'YData');
		theTag = get(h, 'Tag');
		if isequal(theTag, theCornerTag)
			theCorners = [theCorners k];
		end
	end
else
	x = thePoints(:, 1); x = x(:).';
	y = thePoints(:, 2); y = y(:).';
	t = thePoints(:, 3); t = t(:).';
	theCorners = find(t);
end

z = x + sqrt(-1)*y;

% Do gridding along splines.  The independent
%  variable is index-number, rather than distance
%  along the curve.  Iterated splines based on the
%  latter are not asymptotically stable.

xtemp = x;
xtemp(end+1) = xtemp(1);
ytemp = y;
ytemp(end+1) = ytemp(1);
ztemp = xtemp + sqrt(-1)*ytemp;
ctemp = theCorners;
ctemp(end+1) = length(xtemp);
% m = theGridSize(1);
% n = theGridSize(2);
m = theGridSize(1)+1;
n = theGridSize(2)+1;
count = [m n m n];

theSpacings = psget(self, 'itsSpacings');
if isempty(theSpacings)
	theSpacings{1} = linspace(0, 1, m);
	theSpacings{2} = linspace(0, 1, n);
	psset(self, 'itsSpacings', theSpacings)   % <== Note.
end

theSpacedEdges = psget(self, 'itsSpacedEdges');
s1 = theSpacings{1};
s2 = theSpacings{2};
if theSpacedEdges(1) == 3
	s1 = 1 - s1(end:-1:1);
end
if theSpacedEdges(2) == 4
	s2 = 1 - s2(end:-1:1);
end

zti = [];

for k = 1:4
	cti = linspace(ctemp(k), ctemp(k+1), count(k));
	if useMexFile
		switch k
		case 1
			s = s1;
		case 2
			s = s2;
		case 3
			s = s1;
			s = 1 - s(end:-1:1);
		case 4
			s = s2;
			s = 1 - s(end:-1:1);
		otherwise
			error([' ## Bad switch index: k = ' int2str(k)])
		end
		cti = ctemp(k) + s * (ctemp(k+1) - ctemp(k));
	end
	
	cti = cti - min(cti);
	cti = cti / max(cti);
	
	ct = ctemp(k):ctemp(k+1);
	ct = ct - min(ct);
	ct = ct / max(ct);
	
	zt = ztemp(ctemp(k):ctemp(k+1));
	
	pptemp = feval(theInterpFcn, ct(:), zt(:), cti, theInterpMethod);
	
	dist = [0; cumsum(abs(diff(pptemp)))];
	dist = dist - min(dist);
	dist = dist / max(dist);
	
	pptemp = feval(theInterpFcn, dist, pptemp, cti, theInterpMethod);
	pptemp(end) = [];

	zti = [zti pptemp(:).'];   % Concession to Unix Matlab.
end

ci = cumsum([1 m-1 n-1 m-1]);

% Roll to activate edges 3 and/or 4.

theSpacedEdges = psget(self, 'itsSpacedEdges');
if isempty(theSpacedEdges)
	theSpacedEdges = [1 2];
	psset(self, 'itsSpacedEdges', theSpacedEdges)
end

theRoll = 0;
switch theSpacedEdges(1)
case 1
	switch theSpacedEdges(2)
	case 2
		theRoll = 0;
	case 4
		theRoll = 3;
	end
case 3
	switch theSpacedEdges(2)
	case 2
		theRoll = 1;
	case 4
		theRoll = 2;
	end
end

isCorner = zeros(size(zti));
isCorner(ci) = 1;
for i = 1:theRoll
	f = find(isCorner);
	f = f(2);
	zti = [zti(f:end) zti(1:f-1)];
	isCorner = [isCorner(f:end) isCorner(1:f-1)];
	temp = m;
	m = n;
	n = temp;
end
ci = find(isCorner);

% Compute the grid.

theFcn = 'rect2grid';
if useMexFile, theFcn = 'mexrect2grid'; end

[w, err] = feval(theFcn, zti, [], ci, [m n]);   % <===

if (0)
	disp([' ## ' mfilename])
	gridcalc_laplacians = err;
	disp([' ## gridcalc_laplacians = ' num2str(err)])
end   % Laplacian error-norm.

% Roll back to original orientation.

for i = 1:theRoll
	w = flipud(w.');
end

theGrids = {real(w), imag(w)};

psset(self, 'itsGrids', theGrids)

if nargout > 0, theResult = self; end
