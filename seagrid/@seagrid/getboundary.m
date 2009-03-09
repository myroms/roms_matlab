function theResult = getboundary(self, theCommand, theMessage)

% seagrid/getboundary -- Boundary construction.
%  getboundary(self, theCommand, theMessage) processes theCommand
%   and theMessage on behalf of self, a "seagrid" object.
%
%  [> CORNERS <] -- Start by clicking on four locations in
%     counter-clockwise sequence to initialize the four
%     corner-points of the boundary.
%
%  [> EDGES <] -- Use "mousedown-and-drag" on any edge to insert
%     and drag a new control-point.
%
%  [> DRAG <] -- Use "mousedown-and-drag" on a control point or
%     spacer-point to move it.
%
%  [> DELETE <] -- Use "alt-click" on a control-point to delete it.
%
%  [> ZOOM <] -- Click in white-space to invoke "zoomsafe" action.
%     See "help zoomsafe".
 
% svn $Id$
%=======================================================================
% Copyright (C) 1999 Dr. Charles R. Denham, ZYDECO.
%  All Rights Reserved.
%   Disclosure without explicit written consent from the
%    copyright owner does not constitute publication.
%=======================================================================
 
% Version of 09-Apr-1999 09:19:32.
% Updated    23-Apr-2001 16:24:11.

persistent theDraggable
persistent didMove

if nargin < 1, help(mfilename), return, end
if nargin < 2, theCommand = 'start'; end
if nargin < 3, theMessage = ''; end

% Initialization.

switch lower(theCommand)
case 'start'
	h = get(gca, 'Children');
	if isempty(h)
		set(gca, 'XLim', [-0.5 0.5], 'YLim', [-0.5 0.5])
	else
		axis tight
	end
	doticks(self)
	theButtonDownFcn = 'psevent corner';
	set(gcf, 'WindowButtonDownFcn', theButtonDownFcn)
	if ismac
		set(gcf, 'Pointer', 'crosshair')
	else
		set(gcf, 'Pointer', 'cross')
	end
%	help_start(self)
	return
case 'clear'
	set(gca, 'ButtonDownFcn', '', 'UserData', [])
	return
end

% Mouse activities.

% disp([' ## ' command ' ' num2str(gcbo)])

curr = get(gca, 'CurrentPoint');
x = curr(1, 1); y = curr(1, 2);
% disp([' ## ' num2str(x) ' ' num2str(y)])
sel = get(gcf, 'SelectionType');

theLineStyle = psget(self, 'itsCornerLineStyle');
theMarker = psget(self, 'itsCornerMarker');
theColor = psget(self, 'itsCornerColor');
theEraseMode = psget(self, 'itsEraseMode');

% Get four corner-points.

switch lower(theCommand)
case 'corner'
	thePoints = get(gca, 'UserData');
	if length(thePoints) < 4
		theTag = 'corner-point';
		hold on
		h = plot(x, y, 'EraseMode', theEraseMode, 'Tag', theTag);
		hold off
		thePoints = [thePoints h];
		switch length(thePoints)
		case 1
			theMarker = '*';
		otherwise
			theMarker = 'o';
		end
		set(h, 'LineStyle', theLineStyle, 'Marker', theMarker, ...
				'Color', theColor, 'EraseMode', theEraseMode)
		set(gca, 'UserData', thePoints)
	end
	if length(thePoints) >= 4
		set(gcf, 'Pointer', 'arrow')
		thePoints = get(gca, 'UserData');
		thePoints = rearrange(thePoints(1:4));
		set(gca, 'ButtonDownFcn', '', 'UserData', [])
		set(gcf, 'WindowButtonDownFcn', '')
		h = findobj('Type', 'line', 'Tag', 'coastline');
		if any(h), set(h, 'ButtonDownFcn', ''), end
		h = findobj('Type', 'line', 'Tag', 'bathymetry');
		if any(h), set(h, 'ButtonDownFcn', ''), end
		psset(self, 'itsPoints', thePoints)
		psset(self, 'itsMaskToolFlag', 0)   % Note.
		doupdate(self, 1)
	end
	return
otherwise
end

% Manipulate points.

switch lower(theCommand)
case 'motion'
	didMove = 1;
	switch get(theDraggable, 'Tag')
	case {'corner-point', 'edge-point', 'pot'}
		set(theDraggable, 'XData', x, 'YData', y)
	case 'spacer'
		switch sel
		case 'alt'
			set(theDraggable, 'XData', x, 'YData', y)
		otherwise
			set(theDraggable, 'XData', x, 'YData', y)
		end
	otherwise
	end
case 'up'
	set(gcf, 'WindowButtonMotionFcn', '', 'WindowButtonUpFcn', '')
	set(gcf, 'Pointer', 'arrow')
	switch get(theDraggable, 'Tag')
	case {'corner-point', 'edge-point', 'pot'}
		set(theDraggable, 'XData', x, 'YData', y)
		theDraggable = [];
		
		psset(self, 'itsMaskToolFlag', 0)
		
		doupdate(self, 1)
	case 'spacer'
		set(theDraggable, 'XData', x, 'YData', y)
		theDraggable = [];
		sel = get(gcf, 'SelectionType');
		switch sel
		case 'normal'
			psset(self, 'itsMaskToolFlag', 0)   % Note.
			doupdate(self, ~ismac)
		case 'extend'
			dospacings(self)
			dospacers(self)
		case 'alt'
		end
	otherwise
	end
	didMove = 0;
case 'down'
	didMove = 0;
	switch get(gco, 'Tag')
	case 'corner-point'   % Drag a corner-point.
		switch sel
		case 'normal'   % Drag it.
			theDraggable = gco;
			getboundary(self, 'dodrag')
		otherwise
			disp([' ## "' sel '-click" not used by corners.'])
		end
	case 'edge'   % Add/drag a new edge-point.
		switch sel
		case 'normal'   % Add/drag it.
			theEraseMode = psget(self, 'itsEraseMode');
			hold on
			h = plot(x, y, 'EraseMode', theEraseMode, 'Tag', 'edge-point');
			hold off
			theDraggable = h;
			thePoints = psget(self, 'itsPoints');
			p = get(gco, 'UserData');
			p(~isnan(p)) = thePoints;
			p(isnan(p)) = h;
			thePoints = p;
			theButtonDownFcn = 'psevent down';
			theLineStyle = 'none';
			theMarker = psget(self, 'itsEdgeMarker');
			theColor = psget(self, 'itsEdgeColor');
			set(h, 'UserData', thePoints, 'ButtonDownFcn', theButtonDownFcn)
			set(h, 'EraseMode', theEraseMode, 'LineStyle', 'none', ...
					'Marker', theMarker, 'Color', theColor)
			psset(self, 'itsPoints', thePoints)
			getboundary(self, 'dodrag')
		otherwise
			disp([' ## "' sel '-click" not used by edges.'])
		end
	case 'edge-point'   % Drag or delete edge-point.
		switch sel
		case 'normal'   % Drag it.
			theDraggable = gco;
			didMove = 0;
			getboundary(self, 'dodrag')
		case 'extend'   % Toggle -- not used.
			theSelected = get(gco, 'Selected');
			switch theSelected
			case 'on'
				theSelected = 'off';
			case 'off'
				theSelected = 'on';
			otherwise
			end
			set(gco, 'Selected', theSelected)
		case 'alt'   % Delete it.
			thePoints = get(gco, 'UserData');
			thePoints(thePoints == gco) = [];
			delete(gco)
			set(thePoints, 'UserData', thePoints)
			psset(self, 'itsPoints', thePoints)
			psset(self, 'itsMaskToolFlag', 0)   % Note.
			doupdate(self, 1)
		otherwise
		end
		
% N.B. We will double or (halve) the local grid density on
%  either side of a spacer by using click (or shift-slick).
%  Option-click can be used to reposition a single grid-line.

	case 'spacer'
		switch sel
		case 'normal'   % Select it.
			h = findobj('Tag', 'pot', 'Selected', 'on');
			set(h, 'Selected', 'off')
			theDraggable = gco;
			set(theDraggable, 'Selected', 'on')
			didMove = 0;
			getboundary(self, 'dodrag')
		case 'extend'   % Select it.
			h = findobj('Tag', 'pot', 'Selected', 'on');
			set(h, 'Selected', 'off')
			theDraggable = gco;
			set(theDraggable, 'Selected', 'on')
			didMove = 0;
			getboundary(self, 'dodrag')
		case 'alt'   % Drag it.
			h = findobj('Tag', 'pot', 'Selected', 'on');
			set(h, 'Selected', 'off')
			theDraggable = gco;
			set(theDraggable, 'Selected', 'on')
			didMove = 0;
			getboundary(self, 'dodrag')
		otherwise
			disp([' ## "' sel '-click" not used by spacers.'])
		end
	case 'pot'
	otherwise
	end
case 'dodrag'   % Drag one point.
	thePointer = 'circle';
	theSelected = 'on';
	theWindowButtonMotionFcn = 'psevent motion';
	theWindowButtonUpFcn = 'psevent up';
	set(gcf, 'WindowButtonMotionFcn', theWindowButtonMotionFcn, ...
				'WindowButtonUpFcn', theWindowButtonUpFcn, ...
				'Pointer', thePointer)
otherwise
end

function theResult = rearrange(thePoints)

% rearrange -- Rearrange points to counter-clockwise.
%  rearrange(thePoints) flips the order of thePoints
%   (handles) if they are not in counter-clockwise
%   sequence.
 
% Copyright (C) 1999 Dr. Charles R. Denham, ZYDECO.
%  All Rights Reserved.
%   Disclosure without explicit written consent from the
%    copyright owner does not constitute publication.
 
% Version of 16-Apr-1999 11:19:36.

if nargin < 1, help(mfilename), return, end

x = zeros(length(thePoints), 1);
y = zeros(size(x));

for i = 1:length(thePoints)
	x(i) = get(thePoints(i), 'XData');
	y(i) = get(thePoints(i), 'YData');
end

x(end+1) = x(1);
y(end+1) = y(1);

i = 1:length(x)-1;
j = i+1;

theArea = 0.5 * sum(x(i).*y(j) - x(j).*y(i));

result = thePoints;
if theArea < 0
	result(2:length(result)) = result(length(result):-1:2);
end

theResult = result;
