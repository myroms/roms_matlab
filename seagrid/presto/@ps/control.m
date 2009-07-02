function theResult = control(self, theStyle, ...
						theNormalizedPosition, thePixelOffset, ...
						varargin)

% ps/control -- Create UIControl for "ps".
%  control(self, 'theStyle', theNormalizedPosition, thePixelOffset)
%   creates the given control and enables it to be repositioned
%   in normalized and pixel coordinates during resizing.  The
%   style may be given as a string, or as an actual control
%   already present in the figure associated with self, a "px"
%   object.  Additional name/value pairs can be appended to the
%   argument list to change the properties of the control.
%
%   Note: scrollbars along the edge of the figure can be specified
%   by their location: 'bottom', 'right', 'top', or 'left'.  The
%   "String" property will be set to one of those.  For other
%   controls, the "String" can be set with a name/value pair.
 
% svn $Id$
%=======================================================================
% Copyright (C) 1999 Dr. Charles R. Denham, ZYDECO.
%  All Rights Reserved.
%   Disclosure without explicit written consent from the
%    copyright owner does not constitute publication.
%=======================================================================
 
% Version of 03-Nov-1999 00:34:10.
% Updated    04-Nov-1999 16:05:09.

if nargout > 0, theResult = []; end
if nargin < 1, help(mfilename), return, end
if nargin < 2, theStyle = 'pushbutton'; end
if nargin < 3
	theNormalizedPosition = [0 0 0 0];
end
if nargin < 4
	thePixelOffset = get(0, 'DefaultUIControlPosition');
end

theTag = '';

switch lower(theStyle)
case {'xscrollbar', 'hscrollbar', 'bottom'}
	theString = 'bottom';
	theTag = 'bottom';
	theStyle = 'slider';
	theNormalizedPosition = [0 0 1 0];
	thePixelOffset = [16 0 -32 16];
case {'yscrollbar', 'vscrollbar', 'right'}
	theString = 'right';
	theTag = 'right';
	theStyle = 'slider';
	theNormalizedPosition = [1 0 0 1];
	thePixelOffset = [-16 16 16 -32];
case {'zscrollbar', 'left'}
	theString = 'left';
	theTag = 'left';
	theStyle = 'slider';
	theNormalizedPosition = [0 0 0 1];
	thePixelOffset = [0 16 16 -32];
case {'cscrollbar', 'top'}
	theString = 'top';
	theTag = 'top';
	theStyle = 'slider';
	theNormalizedPosition = [0 1 1 0];
	thePixelOffset = [16 -16 -32 16];
end

theLayout = [theNormalizedPosition; thePixelOffset];
u.itsLayout = theLayout;

if ishandle(theStyle)
	theControl = theStyle;
else
	theControl = uicontrol('Style', theStyle);
	if ~isempty(varargin)
		set(theControl, varargin{:})
	end
	set(theControl, 'UserData', u, 'Tag', theTag)
end

doresize(self)

if nargout > 0
	theResult = theControl;
end
