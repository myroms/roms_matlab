function theResult = dograticule(self, isVisible)

% dograticule -- Redraw geographic graticule in SeaGrid.
%  dograticule(self) redraws the map graticule in the
%   SeaGrid window represented by self, a "seagrid" object.
%  domgraticule(self, 'isVisible') shows the map grid
%   if 'isVisible' = 'on'.
 
% svn $Id$
%=======================================================================
% Copyright (C) 2000 Dr. Charles R. Denham, ZYDECO.
%  All Rights Reserved.
%   Disclosure without explicit written consent from the
%    copyright owner does not constitute publication.
%=======================================================================
 
% Version of 06-Jan-2000 14:14:03.
% Updated    07-Jan-2000 07:12:45.

theFigure = ps(self);

f = findobj(theFigure, 'Type', 'line', 'Tag', 'graticule');
if any(f), delete(f), end

if nargin < 2
	isVisible = psget(self, 'itsGraticule');
end

switch isVisible
case 'on'

	x = get(gca, 'XTick');
	y = get(gca, 'YTick');
	xlim = get(gca, 'XLim');
	ylim = get(gca, 'YLim');
	
	x = [xlim(1) x xlim(2)];
	y = [ylim(1) y ylim(2)];
	
	[lon, ignore] = sg_xy2ll(x, min(ylim)*ones(size(x)));
	[ignore, lat] = sg_xy2ll(min(xlim)*ones(size(y)), y);
	
	lon = interp1(linspace(0, 1, length(lon)), lon, ...
					linspace(0, 1, 2*length(lon)-1));
	
	lat = interp1(linspace(0, 1, length(lat)), lat, ...
					linspace(0, 1, 2*length(lat)-1));
	
	[lon, lat] = meshgrid(lon, lat);
	
	[x, y] = sg_ll2xy(lon, lat);
	
	hold on
	for i = 1:2
		plot(x(:, 3:2:end-2), y(:, 3:2:end-2), 'b:', ...
			'EraseMode', 'xor', 'Tag', 'graticule')
		if i == 1, x = x.'; y = y.'; end
	end
	hold off
	
case 'off'
otherwise
end
