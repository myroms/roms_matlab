function theResult = doscribble(self, theEvent, varargin)

% presto/doscribble -- Process "presto" mouse events.
%  doscribble(self, 'theEvent') handles mouse events
%   on behalf of self, a "presto" object.  The mouse
%   track is temporarily scribbled in the window.

% svn $Id$
%=======================================================================
% Copyright (C) 1999 Dr. Charles R. Denham, ZYDECO.
%  All Rights Reserved.
%   Disclosure without explicit written consent from the
%    copyright owner does not constitute publication.
%=======================================================================
 
% Version of 05-Nov-1999 08:43:30.
% Updated    14-Dec-1999 16:46:14.

persistent OLD_NAME
persistent OLD_NUMBER_TITLE
persistent OLD_POINTER
persistent SELECTION_TYPE

switch lower(theEvent)

% Mouse play.
	
case 'windowbuttonmotionfcn'
	pt = get(gca, 'CurrentPoint');
	pt = mean(pt);
	NEW_NAME = ['[' num2str(pt(1)) ', ' num2str(pt(2)) ', ' num2str(pt(3)) ']'];
	set(gcbf, 'Name', NEW_NAME)
	h = findobj(gcbf, 'Type', 'line', 'Tag', 'presto-scribble');
	if any(h)
		h = h(1);
		x = get(h, 'XData');
		y = get(h, 'YData');
		x(end+1) = pt(1);
		y(end+1) = pt(2);
		while length(x) > 10
			x(1) = [];
			y(1) = [];
		end
		set(h, 'XData', x, 'YData', y)
	end
case 'windowbuttondownfcn'
	pt = get(gca, 'CurrentPoint');
	pt = mean(pt);
	SELECTION_TYPE = get(gcbf, 'SelectionType');
	NEW_NAME = ['[' num2str(pt(1)) ', ' num2str(pt(2)) ', ' num2str(pt(3)) ']'];
	NEW_NUMBER_TITLE = 'off';
	NEW_POINTER = 'circle';
	OLD_NAME = get(gcbf, 'Name');
	OLD_NUMBER_TITLE = get(gcbf, 'NumberTitle');
	OLD_POINTER = get(gcbf, 'Pointer');
	set(gcbf, ...
			'WindowButtonDownFcn', '', ...
			'WindowButtonMotionFcn', ['event WindowButtonMotionFcn'], ...
			'WindowButtonUpFcn', ['event WindowButtonUpFcn'], ...
			'Name', NEW_NAME, 'NumberTitle', NEW_NUMBER_TITLE, ...
			'Pointer', NEW_POINTER);
	h = line(pt(1), pt(2), 'EraseMode', 'xor', ...
							'Color', [0 0 0], ...
							'LineWidth', 3.0, ...
							'LineStyle', '--', ...
							'Tag', 'presto-scribble');
case 'windowbuttonupfcn'
	set(gcbf, ...
			'WindowButtonMotionFcn', '', ...
			'WindowButtonUpFcn', '', ...
			'WindowButtonDownFcn', ['event WindowButtonDownFcn'], ...
			'Name', OLD_NAME, 'NumberTitle', OLD_NUMBER_TITLE, ...
			'Pointer', OLD_POINTER);
	OLD_NAME = [];
	OLD_POINTER = [];
	h = findobj(gcbf, 'Type', 'line', 'Tag', 'presto-scribble');
	if ishandle(h), delete(h), end
end

if nargout > 0
	theResult = self;
end
