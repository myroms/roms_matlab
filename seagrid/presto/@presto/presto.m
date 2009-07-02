function theResult = presto

% presto/presto -- Constructor for "presto" class.
%  presto(theHandle) attaches a new "presto" object
%   to theHandle, which defaults to a new figure.
%   However, if theHandle already has a "presto"
%   object, that object is simply returned and
%   no new object is created.
%  presto(theHandle) returns the "presto" object
%   already associated with theHandle.
%  presto(thePresto) updates the "presto" object
%   by writing thePresto to its associated handle.
%
% Synopsis of Presto
%
%    Presto is a class for managing events generated
% by graphical objects.  Each "presto" object is attached
% to a Matlab handle.  The properties of the object and
% of its handle are manipulated in the same manner, using
% "dot" syntax.  This allows the user to define arbitrary
% properties in the object, as well as to manipulate the
% graphical properties of the associated handle easily.
% Subscripting can be nested, which makes the "presto"
% syntax more versatile than that used for the Matlab
% "set" and "get" commands.
%
%    Presto preserves any existing "UserData" and allows
% it to be manipulated via self.UserData.  Special handling
% takes place where "UserData" are concerned, because the
% "presto" object is stored there.  Any existing "UserData"
% is preserved and can even be manipulated as "self.UserData".
% When the "presto" object is disconnected via its "detach"
% method, the "UserData" are copied back into the handle.
%
%    Graphical events are expected to call "event" with
% the name of the particular callback.  The "presto/enable"
% method is useful for automatically enabling all possible
% callbacks associated with the handle.  In the case of a
% figure, for example, this would include menus, uicontrols,
% axes, and all other interactive graphical types.
%
%    The users "event" method is responsible for performing
% all processing associated with a particular event.  The
% "presto" base-class pays attention only to a "quit" event,
% due to a "Quit" menu-item, or to clicking in the go-away box.
% Other events that are passed to the base-class are simply
% displayed without further action.
%
%    Uses can invoke the "inherit" command to pass events
% to the base-class; its syntax is similar to "feval".
%
%    Most activities can be handled by a single "presto"
% object attached to a figure.  Some situations may warrant
% other arrangements, however.  For example, the resizing
% of a window can generate disarray amongst the controls,
% for which an automated "layout" scheme may be needed.
% It would make sense to create a "layout" object for
% each control, whose sole method would be "resize",
% to be called by the figure whenever it noticed a
% "ResizeFcn" event.
%
%    Presto reacts to the following events,
% which the derived class should trap as needed:
%
% new, open, close, save, saveas,
% revert, pagesetup, print, quit,
% undo, cut, copy, paste, clear
% windowbuttondownfcn, windowbuttonmotionfcn, windowbuttonupfcn,
% resizefcn, buttondownfcn, callback
% closerequestfcn, deletefcn

% svn $Id$
%=======================================================================
% Copyright (C) 1999 Dr. Charles R. Denham, ZYDECO.
%  All Rights Reserved.
%   Disclosure without explicit written consent from the
%    copyright owner does not constitute publication.
%=======================================================================
 
% Version of 27-Oct-1999 23:08:50.
% Updated    05-Nov-1999 20:52:48.

theFigure = figure('Name', 'Presto Demo', 'MenuBar', 'none');
theStruct.ignore = [];
self = class(theStruct, 'presto', ps(theFigure));
psbind(self)

theMenus = { ...
	'<Presto Demo>', ...
	'>File', '>>MenuBar', '->Quit', ...
	'>Help', '>>Controls', '>>Events', ...
	'>>Handlers', '>>Menus', '>>Subclassing', ...
	'->Version', ...
};

menu(self, theMenus)

h = control(self, 'bottom');
set(h, 'Min', 0, 'Max', 1, 'Value', 0.5);
h = control(self, 'right');
set(h, 'Min', 0, 'Max', 1, 'Value', 0.5);
h = control(self, 'left');
set(h, 'Min', 0, 'Max', 1, 'Value', 0.5);
h = control(self, 'top');
set(h, 'Min', 0, 'Max', 1, 'Value', 0.5);

enable(self);

theEventHandlers = { ...
	'WindowButtonMotionFcn', 'doscribble', ...
	'WindowButtonDownFcn', 'doscribble', ...
	'WindowButtonUpFcn', 'doscribble', ...
	'ResizeFcn', 'doresize', ...
	'bottom', 'doroll', ...
	'right', 'doroll', ...
	'left', 'doroll', ...
	'top', 'doroll', ...
	'menubar', 'domenubar', ...
	'controls', 'help_controls', ...
	'events', 'help_events', ...
	'handlers', 'help_handlers', ...
	'menus', 'help_menus', ...
	'subclassing', 'help_subclassing', ...
	'version', 'help_version', ...
	'quit', 'doquit', ...
	'CloseRequestFcn', 'doquit', ...
};

handler(self, theEventHandlers{:});

n = 36;
n = n + rem(n, 4);
[x, y, z] = sphere(n);
x([1 end], :) = []; y([1 end], :) = []; z([1 end], :) = [];
c = 1 + 0.1*rand(size(x));
c = zeros(size(x));
c(:, 1:n/4) = 1;
c(:, n/2+1:n/2+n/4) = 1;
c(:, end) = c(:, 1);
surface(x, y, z, c)
xlabel x, ylabel y, zlabel z
set(gca, 'Color', get(gcf, 'Color'))
set(gca, 'Visible', 'off')
view(0, 0)
axis equal
disp(' ## Click and move mouse in "presto" window.')
if nargout > 0
	theResult = self;
else
	assignin('caller', 'ans', self)
end

% Get the "presto" object from an existing handle.

if nargin > 0
	if ishandle(theHandle)
		self = ps(theHandle);
		if isps
			if nargout > 0
				theResult = self;
			else
				assignin('caller', 'ans', self)
				disp(self)
			end
			return
		end
	elseif isps(theHandle, 'ps')
		self = theHandle;
		if nargout > 0
			theResult = self;
		else
			assignin('caller', 'ans', self)
			disp(self)
		end
		return
	end
end

if nargout > 0
	theResult = self;
else
	assignin('caller', 'ans', self)
end
