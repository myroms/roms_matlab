function theResult = doticks(self, autoBounds)

% seagrid/doticks -- Geographic tick-marks.
%  doticks(self) updates the geographic tickmarks
%   on behalf of self, a seagrid object.
%  doticks(self, autoBounds) causes the geographic
%   bounds to be set automatically if "autoBounds"
%   is logically TRUE, even if manual bounds are
%   available.
 
% svn $Id$
%=======================================================================
% Copyright (C) 1999 Dr. Charles R. Denham, ZYDECO.
%  All Rights Reserved.
%   Disclosure without explicit written consent from the
%    copyright owner does not constitute publication.
%=======================================================================
 
% Version of 26-Jul-1999 10:45:47.
% Updated    11-Feb-2000 16:28:32.

if nargout > 0, theResult = self; end
if nargin < 1, help(mfilename), return, end
if nargin < 2, autoBounds = 0; end

if verbose(self)
	disp([' ## ' mfilename])
end

theProjection = psget(self, 'itsProjection');
theMapUnits = psget(self, 'itsMapUnits');

theFigure = ps(self);
setsafe(0, 'CurrentFigure', theFigure)

% Project plot limits to geographic coordinates.

theAxes = findobj('Type', 'axes', 'Tag', class(self));
if isempty(theAxes)
	theAxes = gca;
	set(theAxes, 'Tag', class(self))
end

setsafe(theFigure, 'CurrentAxes', theAxes)

setsafe(theAxes, 'DataAspectRatio', [1 1 1])

xlim = get(theAxes, 'XLim');
ylim = get(theAxes, 'YLim');
axpos = get(theAxes, 'Position');

axpos(1) = axpos(1) + 10*axpos(3);   % Push out of window.

% Note: bounds that are manually set override any zooming
%  that could have taken place.  We need to figure out how
%  to have them both.

theLongitudeBounds = psget(self, 'itsLongitudeBounds');
theLatitudeBounds = psget(self, 'itsLatitudeBounds');
	
if length(theLongitudeBounds) == 2 & length(theLatitudeBounds) == 2
	if ~autoBounds
		lo = theLongitudeBounds;
		la = min(theLatitudeBounds) * [1 1];
		[xlim, ignore] = sg_ll2xy(lo, la);
		lo = min(theLongitudeBounds) * [1 1];
		la = theLatitudeBounds;
		[ignore, ylim] = sg_ll2xy(lo, la);
		setsafe(theAxes, 'XLim', xlim, 'YLim', ylim)
	else
		axis normal
	end
elseif length(theLongitudeBounds) ~= length(theLatitudeBounds)
	disp(' ## Longitude and latitude bounds require 2 elements each.')
	axis normal
else
if verbose(self), hello doticks line 79 axis normal disabled, end
%	axis normal
end

switch lower(theMapUnits)

case {'degrees', 'deg'}
	[lolim, lalim] = sg_xy2ll(xlim, ylim);

	if ~isunix | 1

% Use a temporary, invisible axis to get tick-marks
%  and labels for the equivalent geographic ranges.
%  This is a poorman's "linticks" scheme.

		ax = axes('Position', axpos, 'Visible', 'off', ...
					'XLim', sort(lolim), 'YLim', sort(lalim));
		
		lotick = get(ax, 'XTick');
		loticklabel = get(ax, 'XTickLabel');
		
		latick = get(ax, 'YTick');
		laticklabel = get(ax, 'YTickLabel');
		
		delete(ax)
	else
		
% Try a "linticks" scheme that avoids creating any
%  additional graphics, even invisible ones.

		lotick = linticks(min(lolim), max(lolim), 10);
		f = find(lotick >= min(lolim) & lotick <= max(lolim));
		lotick = lotick(f);
		loticklabel = num2str(lotick(:));   % Right-justified.
		[m, n] = size(loticklabel);
		for i = 1:m
			lab = loticklabel(i, :);
			k = floor(sum(lab == ' ')/2);
			lab = [lab(k+1:end) lab(1:k)];   % Centered.
			loticklabel(i, :) = lab;
		end
		
		latick = linticks(min(lalim), max(lalim), 10);
		f = find(latick >= min(lalim) & latick <= max(lalim));
		latick = latick(f);
		laticklabel = num2str(latick(:));   % Right-justified.
	end

% Place the geographic ticks and labels
%  on the projected seagrid map.

	lo = lotick;
	la = latick(1) + zeros(size(lotick));
	[xtick, ignore] = sg_ll2xy(lo, la);
	
	lo = lotick(1) + zeros(size(latick));
	la = latick;
	[ignore, ytick] = sg_ll2xy(lo, la);
	
	setsafe(theAxes, 'XTick', xtick, 'XTickLabel', loticklabel, ...
						'YTick', ytick, 'YTickLabel', laticklabel)

	xlab = get(theAxes, 'XLabel');
	ylab = get(theAxes, 'YLabel');
	setsafe(xlab, 'String', 'Longitude')
	setsafe(ylab, 'String', 'Latitude')
	
case {'kilometers', 'km'}
	EARTH_RADIUS = 6371.315;   % kilometers.
	xlim = xlim * EARTH_RADIUS;
	ylim = ylim * EARTH_RADIUS;
	ax = axes('Position', axpos, 'Visible', 'off', ...
				'XLim', sort(xlim), 'YLim', sort(ylim));
	xtick = get(ax, 'XTick') / EARTH_RADIUS;
	xticklabel = get(ax, 'XTickLabel');
	ytick = get(ax, 'YTick') / EARTH_RADIUS;
	yticklabel = get(ax, 'YTickLabel');
	delete(ax)
	setsafe(theAxes, 'XTick', xtick, 'XTickLabel', xticklabel, ...
						'YTick', ytick, 'YTickLabel', yticklabel)

	xlab = get(theAxes, 'XLabel');
	ylab = get(theAxes, 'YLabel');
	setsafe(xlab, 'String', 'X (km)')
	setsafe(ylab, 'String', 'Y (km)')
	
case 'projected'
	xlim = xlim;
	ylim = ylim;
	ax = axes('Position', axpos, 'Visible', 'off', ...
				'XLim', sort(xlim), 'YLim', sort(ylim));
	xtick = get(ax, 'XTick');
	xticklabel = get(ax, 'XTickLabel');
	ytick = get(ax, 'YTick');
	yticklabel = get(ax, 'YTickLabel');
	delete(ax)
	setsafe(theAxes, 'XTick', xtick, 'XTickLabel', xticklabel, ...
						'YTick', ytick, 'YTickLabel', yticklabel)

	xlab = get(theAxes, 'XLabel');
	ylab = get(theAxes, 'YLabel');
	setsafe(xlab, 'String', 'X')
	setsafe(ylab, 'String', 'Y')
otherwise
end

dograticule(self)

setsafe(gca, 'DataAspectRatio', [1 1 1])

% Would it be helpful to do "drawnow" here?  No.
% drawnow

if nargout > 0, theResult = self; end


% ---------- linticks ---------- %

function theResult = linticks(theMin, theMax, theTickCount)

% linticks -- Tick positions within a range.
%  linticks(theMin, theMax, theTickCount) returns a vector
%   of tick positions that span the interval from theMin
%   to theMax, with theTickCount or fewer elements.
%   Based on the "goodscales" scheme used by the Matlab
%   "plotyy" routine.
 
% Copyright (C) 2000 Dr. Charles R. Denham, ZYDECO.
%  All Rights Reserved.
%   Disclosure without explicit written consent from the
%    copyright owner does not constitute publication.
 
% Version of 08-Feb-2000 21:27:14.
% Updated    08-Feb-2000 21:27:14.

if nargout > 0, theResult = []; end

if nargin < 1, help(mfilename), return, end

if ischar(theMin), theMin = eval(theMin); end
if ischar(theMax), theMax = eval(theMax); end
if ischar(theTickCount), theTickCount = eval(theTickCount); end

[low, high, ticks] = goodscales(theMin, theMax);

f = find(ticks <= theTickCount);
if any(f)
	f = f(1);
	result = linspace(low(f), high(f), ticks(f));
	result(abs(result) < sqrt(eps)) = 0;
end

if nargout > 0
	theResult = result;
else
	assignin('caller', 'ans', result)
	disp(result)
end


% ---------- goodscales ---------- %


function [low, high, ticks] = goodscales(xmin, xmax)

% GOODSCALES -- Returns parameters for "good" scales.
%  [LOW, HIGH, TICKS] = GOODSCALES(XMIN, XMAX) returns lower and upper
%   axis limits (LOW and HIGH) that span the interval (XMIN, XMAX) 
%   with "nice" tick spacing.  The number of major axis ticks is 
%   returned in TICKS.

% ** Liberated from "plotyy.m".

% Modifications:
% Copyright (C) 2000 Dr. Charles R. Denham, ZYDECO.
%  All Rights Reserved.
%   Disclosure without explicit written consent from the
%    copyright owner does not constitute publication.
 
% Version of 08-Feb-2000 20:50:36.
% Updated    08-Feb-2000 21:16:45.

if nargin < 1, help(mfilename), return, end

if ischar(xmin), xmin = eval(xmin); end
if ischar(xmax), xmax = eval(xmax); end

BestDelta = [ 0.1 0.2 0.5 1 2 5 10 20 50 ];

xmin = min(xmin(:));
xmax = max(xmax(:));

if xmin == xmax
	lo = xmin;
	high = xmax + 1;
	ticks = 1;
else
	Xdelta = xmax-xmin;
	delta = 10.^(round(log10(Xdelta) - 1)) * BestDelta;
	high = delta .* ceil(xmax ./ delta);
	lo = delta .* floor(xmin ./ delta);
	ticks = round((high - lo) ./ delta) + 1;
end

if nargout > 0
	low = lo;
elseif nargout == 1
	low = [lo; high; ticks];
else
	disp([lo; high; ticks])
end

