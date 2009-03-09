function theResult = doedges(self)

% doedges -- Draw the edges.
%  doedges(self) draws the interactive boundary
%   associated with self, a "seagrid" object.
 
% svn $Id$
%=======================================================================
% Copyright (C) 1999 Dr. Charles R. Denham, ZYDECO.
%  All Rights Reserved.
%   Disclosure without explicit written consent from the
%    copyright owner does not constitute publication.
%=======================================================================
 
% Version of 09-Apr-1999 17:03:02.
% Updated    21-Sep-2000 15:02:41.

% N.B. -- We eliminate screen flashing by applying the
%  'xor' EraseMode to each drawing motion and using
%  "plot", rather than "line".

% N.B. -- This routine can be made more efficient by
%  rearranging, rather than eliminating, handles.

if nargout > 0, theResult = self; end

% Initialize.

thePoints = psget(self, 'itsPoints');

if length(thePoints)  < 4, return, end

x = zeros(size(thePoints));
y = zeros(size(thePoints));
theTags = cell(size(thePoints));
theEdges = zeros(size(thePoints));   % Edge handles.
theCorners = [];   % Corner indices.

theEraseMode = psget(self, 'itsEraseMode');
theCornerTag = psget(self, 'itsCornerTag');
theEdgePointTag = psget(self, 'itsEdgePointTag');

% Get the xy data and corner-point indices.

if min(size(thePoints)) == 1
	for k = 1:length(thePoints)
		x(k) = get(thePoints(k), 'XData');
		y(k) = get(thePoints(k), 'YData');
		theTags{k} = get(thePoints(k), 'Tag');
		if isequal(theTags{k}, theCornerTag)
			theCorners = [theCorners k];   % #1 is always a theCorners-point.
		end
	end
else
	x = thePoints(:, 1); x = x(:).';
	y = thePoints(:, 2); y = y(:).';
	t = thePoints(:, 3); t = t(:).';
	theCorners = find(t);
	theTags = cell(size(t));
	for k = 1:length(t)
		switch t(k)
		case 0
			theTags{k} = theEdgePointTag;
		otherwise
			theTags{k} = theCornerTag;
		end
	end
	theEdges = zeros(size(x));   % Edge handles.
end

% Compute splines.
%  Our independent variable is the "index" number, but
%  that leads to seemingly, odd behavior when adding
%  a point.  A better variable would be distance along
%  the edge or curve.  We wonder whether there is a
%  "root" spline, similar to the "root" that attends
%  iterative median-filtering.   Must investigate.

xtemp = x;
xtemp(end+1) = xtemp(1);
ytemp = y;
ytemp(end+1) = ytemp(1);
ztemp = xtemp + sqrt(-1) * ytemp;
ctemp = theCorners;
ctemp(end+1) = length(xtemp);

theEndSlopeFlag = psget(self, 'itsEndSlopeFlag');

for i = 1:4
	ct = ctemp(i):ctemp(i+1);
	zt = ztemp(ctemp(i):ctemp(i+1));
	if theEndSlopeFlag
		s = diff(zt) ./ diff(ct);
		zt = [s(1) zt s(end)];   % End-slopes; see "help spline".
	end
	pp{i} = spline(ct, zt);
end

% Delete old edges.

oldEdges = findobj('Type', 'line', 'Tag', 'edge');
delete(oldEdges)

% Draw new edges.

ppindex = 0;
for k = 1:length(xtemp)-1
	if any(k == ctemp)
		ppindex = ppindex+1;
	end
	xx = linspace(k, k+1, 11);
	zz = ppval(pp{ppindex}, xx);
	hold on
	theEdges(k) = plot(real(zz), imag(zz), 'r-', 'EraseMode', theEraseMode);
	hold off
end

theTag = psget(self, 'itsEdgeTag');
theLineStyle = psget(self, 'itsEdgeLineStyle');

theLineWidth = get(0, 'DefaultLineLineWidth');
theLineWidth = 2.5;

theButtonDownFcn = 'psevent down';

thePointFlag = psget(self, 'itsPointFlag');
if ~thePointFlag, theButtonDownFcn = ''; end

set(theEdges, 'Tag', theTag, 'LineStyle', theLineStyle, ...
				'LineWidth', theLineWidth, ...
				'ButtonDownFcn', theButtonDownFcn)

psset(self, 'itsEdges', theEdges)

% Stash a NaN place-holder into each edge, that will
%  indicate the position of the new point, should
%  that edge be clicked upon.

for k = 1:length(thePoints)
	p = zeros(1, length(thePoints)+1);
	p(k+1) = NaN;
	set(theEdges(k), 'UserData', p)
end

if nargout > 0, theResult = self; end
