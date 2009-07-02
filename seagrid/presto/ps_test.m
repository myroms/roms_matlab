function theResult = ps_test

% ps_test -- Test procedure for "ps" class.
%  ps_test (no argument) creates a window
%   and returns the controlling "ps" object.
%   Note: Not all of the mouse actions have
%   corresponding handlers in this procedure.
 
% svn $Id$
%=======================================================================
% Copyright (C) 1999 Dr. Charles R. Denham, ZYDECO.
%  All Rights Reserved.
%   Disclosure without explicit written consent from the
%    copyright owner does not constitute publication.
%=======================================================================
 
% Version of 28-Oct-1999 09:27:50.
% Updated    14-Dec-1999 16:50:42.

setdef(mfilename)

p = ps(figure('Name', mfilename, 'Units', 'normalized'));
figure(handle(p))

% Handlers.

theEventHandlers = { ...
	'windowbuttonmotionfcn', 'doscribble', ...
	'windowbuttondownfcn', 'doscribble', ...
	'windowbuttonupfcn', 'doscribble', ...
	'resizefcn', 'doresize', ...
	'closerequestfcn', 'doquit',...
	...
	'ps', 'dops', ...
	'about', 'doabout', ...
	'file', 'dofile', ...
	'new', 'donew', ...
	'open', 'doopen', ...
	'close', 'doclose', ...
	'save', 'dosave', ...
	'saveas', 'dosaveas', ...
	'revert', 'dorevert', ...
	'pagesetup', 'dopagesetup', ...
	'print', 'doprint', ...
	'quit', 'doquit', ...
	...
	'edit', 'doedit', ...
	'undo', 'doundo', ...
	'cut', 'docut', ...
	'copy', 'docopy', ...
	'paste', 'dopaste', ...
	'clear', 'doclear', ...
	'home', 'dohome', ...
	'menubar', 'domenubar', ...
	'view', 'doview', ...
	...
	'controls', 'help_controls', ...
	'events', 'help_events', ...
	'handlers', 'help_handlers', ...
	'menus', 'help_menus', ...
	'subclassing', 'help_subclassing', ...
	'version', 'help_version', ...
	...
	'bottom', 'doscroll', ...
	'right', 'doscroll', ...
	'left', 'doscroll', ...
	'top', 'doscroll', ...
};

p = handler(p, theEventHandlers{:});

% Menus.

thePSMenu = {
	'<PS>'
};

theAboutMenu = {
	'About...'
};

theFileMenu = {
	'File'
	'>New...'
	'>Open...'
	'-Close'
	'>Save'
	'>Save As...'
	'>Revert'
	'-Page Setup...'
	'>Print...'
	'-Quit'
};

theEditMenu = {
	'Edit'
	'>Undo'
	'-Cut'
	'>Copy'
	'>Paste'
	'>Clear'
	{'-Home', 'H', 'off'}
	{'>Menubar', 'M', 'on'}
	{'>View', '=', 'off'}
};

theHelpMenu = {
	'Help'
	'>Controls...'
	'>Events...'
	'>Handlers...'
	'>Menus...'
	'>Subclassing...'
	'>Version...'
};

theMenus = [thePSMenu; theAboutMenu; ...
				theFileMenu; theEditMenu; theHelpMenu];
for i = 2:length(theMenus)
	if iscell(theMenus{i})
		theMenus{i}{1} = ['>' theMenus{i}{1}];
	else
		theMenus{i} = ['>' theMenus{i}];
	end
end

menu(p, theMenus)

patch
surface
text(1, .5, 'Hello World! This is PS.')
x = (0.1:0.02:0.9);
y = rand(size(x)) + 1.5;
line(x, y)

h(1) = control(p, 'right');
h(2) = control(p, 'bottom');
h(3) = control(p, 'left');
h(4) = control(p, 'top');
set(h, 'Min', 0, 'Max', 1, 'Value', 0.5)

h = control(p, 'radiobutton', [0.5 0.5 0 0], [-30 -20 60 20]);
set(h, 'String', 'Okay')

p = enable(p);

p = ps(gcf);

p.Color = [0 100 100]/100;
p.MenuBar = 'none';

if nargout > 0
	theResult = p;
else
	assignin('caller', 'ans', p)
end
