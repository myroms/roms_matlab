function theResult = grid2mask(x, y, z, c, mask)

% grid2mask -- Create an interactive mask.
%  grid2mask('demo') demonstrates itself via "grid2mask(10)".
%  grid2mask(N) demonstrates itself with an N-by-N grid.
%
%  grid2mask(x, y, z, c, mask) creates a masked surface from the
%   (x, y, z, c, mask) grid, whose cells can be toggled by mouse
%   clicks.  Clicking button #1 on a grid cell causes its color
%   to toggle to/from transparency via NaN.  Button #2 toggles
%   the value with the defined replacement-depth (see below).
%   Button #3 captures the clicked cell, making that the
%   replacement depth for subsequent clicks with Button #2.
%   The default depth is 0.  The default color array c is z itself.
%   A handle is returned, whose "UserData" contains the handles
%   of surfaces that can be used to re-assemble the color array.
%
%  grid2mask('mask', mask_array) sets/gets the color-array of the mask.
%   If logical, the TRUE values are treated as NaNs.
%
%  grid2mask('depth', depth_array) sets/gets the depth-array.
%
%  grid2mask('replacement', rep_value) sets/gets the replacement
%   depth to be used for subsequent clicks with Button #2,
%   whenever the grid-cell is not masked by NaN.
%
%  grid2mask('edit') captures the edit-field depth replacement.
%
%  grid2mask('done') deletes the grid2mask objects.

% svn $Id$
%=======================================================================
% Copyright (C) 1999 Dr. Charles R. Denham, ZYDECO.
%  All Rights Reserved.
%   Disclosure without explicit written consent from the
%    copyright owner does not constitute publication.
%=======================================================================
 
% Version of 07-Oct-1999 15:33:50.
% Updated    23-Apr-2001 10:59:01.

if nargin < 1, help(mfilename), return, end

if nargin == 1 && isequal(x, 'demo')
	x = 10;
end
if ischar(x) && any(x(1) == '0123456789')
	x = eval(x);
end

% Demonstration.

if nargin < 2 && ~ischar(x)
	f = findobj(gcf, 'Tag', mfilename);
	delete(f)
	nx = x;
	if ischar(nx), nx = eval(nx); end
	ny = nx;
	[x, y] = meshgrid(0:nx, 0:ny);
	z = rand(size(x));
	h = feval(mfilename, x, y, z);
	set(gcf, 'Name', mfilename)
	zoomsafe
	if nargout > 0
		theResult = h;
	else
		assignin('caller', 'ans', h)
	end
	return
end

REPLACEMENT_DEPTH = 0;

MASTER_HANDLE = findobj(gca, 'Type', 'text', 'Tag', mfilename);
if any(MASTER_HANDLE) && ishandle(MASTER_HANDLE(1))
	u = get(MASTER_HANDLE, 'UserData');
	STRIP_HANDLES = u.strips;
	EDIT_HANDLE = u.edit;
	REPLACEMENT_DEPTH = u.depth;
	REPLACEMENT_COLOR = u.color;
end

% Get the edit-field replacement-value.

if isequal(lower(x), 'edit') && nargin == 1
	REPLACEMENT_DEPTH = eval(get(EDIT_HANDLE, 'String'));
	set(EDIT_HANDLE, 'String', num2str(REPLACEMENT_DEPTH));
	u.depth = REPLACEMENT_DEPTH;
	set(MASTER_HANDLE, 'UserData', u)
	return
end

% Set/Get the replacement value.

if isequal(lower(x), 'replacement')
	if nargin > 1
		REPLACEMENT_DEPTH = y(1);
		set(EDIT_HANDLE, 'String', num2str(REPLACEMENT_DEPTH))
		u.depth = REPLACEMENT_DEPTH;
		set(MASTER_HANDLE, 'UserData', u)
	else
		result = REPLACEMENT_DEPTH;
		if nargout > 0
			theResult = result;
		else
			assignin('caller', 'ans', result)
		end
	end
	return
end

% Set the mask.

if isequal(lower(x), 'mask') && nargin > 1
	if nargin < 3 || isempty(y), y = 0; end
	mask = ~~y;
	if length(mask) == 1
		mask = ~~(zeros(size(y)) + mask);
	end
	x = get(u.strips(1), 'XData');
	c = zeros(length(u.strips)+1, length(x));
	for i = 1:length(u.strips)
		c(i:i+1, :) = get(u.strips(i), 'UserData');
	end
	if all(mask(:))
		c(:) = NaN;
	else
		c(mask) = NaN;   % Only does c(1, 1) if all ones.
	end
	for i = 1:length(u.strips)
		set(u.strips(i), 'CData', c(i:i+1, :));
	end
	if nargout > 0
		theResult = isnan(c);
	else
		assignin('caller', 'ans', isnan(c))
	end
	return
end

% Get the mask.

if isequal(lower(x), 'mask') && nargin < 2
	x = get(u.strips(1), 'XData');
	c = zeros(length(u.strips)+1, length(x));
	for i = 1:length(u.strips)
		c(i:i+1, :) = get(u.strips(i), 'CData');
	end
	if nargout > 0
		theResult = isnan(c);
	else
		assignin('caller', 'ans', isnan(c))
	end
	return
end

% Set the depth-array.

if isequal(lower(x), 'depth') && nargin > 1
	return
end

% Get the depth-array.

if isequal(lower(x), 'depth') && nargin < 2
	x = get(u.strips(1), 'XData');
	z = zeros(length(u.strips)+1, length(x));
	for i = 1:length(u.strips)
		z(i:i+1, :) = get(u.strips(i), 'ZData');
	end
	if nargout > 0
		theResult = z;
	else
		assignin('caller', 'ans', z)
	end
	return
end

if isequal(lower(x), 'depth') && nargin < 2
	return
end

% Done.

if nargin == 1 && isequal(lower(x), 'done')
	x = get(u.strips(1), 'XData');
	c = zeros(length(u.strips)+1, length(x));
	for i = 1:length(u.strips)
		c(i:i+1, :) = get(u.strips(i), 'CData');
	end
	delete(STRIP_HANDLES)
	delete(EDIT_HANDLE)
	delete(MASTER_HANDLE)
	if nargout > 0
		theResult = isnan(c);
	else
		assignin('caller', 'ans', isnan(c))
	end
	return
end

% Handle mouse click.

if isequal(lower(x), 'buttondownfcn') && nargin == 1
	selection = get(gcbf, 'SelectionType');
	switch selection
	case {'normal', 'extend', 'alt'}
		h = gcbo;
		x = get(h, 'XData');
		y = get(h, 'YData');
		z = get(h, 'ZData');
		c = get(h, 'CData');
		cx = sum(x); cx = 0.25 * (cx(1:end-1) + cx(2:end));
		cy = sum(y); cy = 0.25 * (cy(1:end-1) + cy(2:end));
		cz = cx + cy*sqrt(-1);
		pt = get(gca, 'CurrentPoint');
		pt = pt(1, 1:2); ptx = pt(1); pty = pt(2);
		ptz = ptx + pty*sqrt(-1);
		dz = abs(cz - ptz);
		f = find(dz == min(dz));
		if any(f)
			f = f(1);
			switch selection
			case 'normal'  % Toggle NaN.
				if isnan(c(1, f))
					c_original = get(h, 'UserData');
					c(1, f) = c_original(1, f);
				else
					c(1, f) = NaN;
				end
				set(h, 'CData', c)
			case 'extend'   % Replace the cell's depth.
				z = get(h, 'ZData');
				z(1, f) = REPLACEMENT_DEPTH;
				set(h, 'ZData', z)
				c = get(h, 'CData');
				c(1, f) = REPLACEMENT_COLOR;
				set(h, 'CData', c)
			case 'alt'   % Capture the cell's depth.
				z = get(h, 'ZData');
				REPLACEMENT_DEPTH = z(1, f);
				c = get(h, 'CData');
				REPLACEMENT_COLOR = c(1, f);
				set(EDIT_HANDLE, 'String', num2str(REPLACEMENT_DEPTH));
				u.depth = REPLACEMENT_DEPTH;
				u.color = REPLACEMENT_COLOR;
				set(MASTER_HANDLE, 'UserData', u)
			end
		end
	otherwise
		switch selection
		case 'open'
			button = 4;
			selection = 'double-click';
		otherwise
			button = NaN;
		end
		disp([' ## ' mfilename ' -- Mouse button not supported: ' int2str(button) ' (' selection ')'])
	end
	return
end

% Setup.

if nargin < 3, z = []; end
if nargin < 4, c = []; end
if nargin < 5, mask = []; end

if isempty(z), z = zeros(size(x)); end

if isempty(c)
	temp = surface(x, y, z, 'Visible', 'off');
	c = get(temp, 'CData');
	delete(temp)
end

if isempty(mask), mask = ~~zeros(size(x)); end

% CAREFUL: When z contains NaN, the corresponding
%  cell is not clickable, unlike when c contains NaN.
%  We should set the corresponding z to zero and
%  adjust c accordingly.

masked_color = c;
masked_color(logical(mask)) = NaN;

theEraseMode = 'xor';

for i = 1:size(x, 1)-1
	STRIP_HANDLES(i) = surface(x(i:i+1, :), y(i:i+1, :), ...
			z(i:i+1, :), masked_color(i:i+1, :), ...
			'UserData', c(i:i+1, :), ...
			'EraseMode', theEraseMode);
end

MASTER_HANDLE = text(x(1), y(1), mfilename, ...
					'Visible', 'off', 'Tag', mfilename);

theBDF = [mfilename ' ButtonDownFcn'];
set(STRIP_HANDLES, 'ButtonDownFcn', theBDF, 'Tag', mfilename)

thePosition = get(0, 'DefaultUIControlPosition').*[0 0 1 1];
thePosition(1) = thePosition(1) + 20;
EDIT_HANDLE = uicontrol('Style', 'edit', 'Tag', mfilename, ...
		'String', num2str(REPLACEMENT_DEPTH), ...
		'Position', thePosition, ...
		'Callback', [mfilename ' edit']);

u.strips = STRIP_HANDLES;
u.edit = EDIT_HANDLE;
u.depth = 0;
u.color = 0;

set(MASTER_HANDLE, 'UserData', u)

if nargout > 0
	theResult = MASTER_HANDLE;
else
	assignin('caller', 'ans', MASTER_HANDLE)
end
