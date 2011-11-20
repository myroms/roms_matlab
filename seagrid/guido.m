function theResult = guido(theStruct, theFigureName, isModal, thePrecision)

% guido -- Get info via Matlab struct.
%  guido('demo') demonstrates itself, using the commands
%   given below.
%  guido(theStruct, 'theFigureName') presents a dialog
%   representing the fieldnames and values of theStruct,
%   a Matlab "struct" which must be compatible with
%   "getinfo" and "setinfo".  If theStruct contains
%   embedded structs, nested dialogs are produced for
%   them as needed.  The empty-matrix [] is returned
%   if the main "Cancel" button is pressed.  Also, use
%   "Cancel" to escape a sub-dialog without saving its
%   most recent changes.  Fields named "help_..." are
%   shown in a separate dialog.  If theFigureName is
%   not given or is empty, the external name of
%   theStruct will be used.
%  guido(theStruct, 'theFigureName', isModal) uses modal
%   dialogs on non-PCWIN machines if "isModal" evaluates
%   to logical(1) TRUE.  The default is to use non-modal
%   dialogs, which allows interactions with other windows,
%   except the Matlab command window.  In either case, the
%   the routine prevents dialogs from being dismissed out
%   of sequence.
%  guido(theStruct, 'theFigureName', isModal, thePrecision)
%   uses thePrecision (default = 4) when displaying numbers
%   in an "edit" control.
%
% Commands used in the demonstration:
%
%   s.help = help(mfilename);
%   s.anEdit = 'some text';
%   s.aNumber = pi;
%   s.aCheckbox = {'checkbox' 0};
%   s.aRadiobutton = {'radiobutton' 0};
%   s.aPopupmenu = {{'red', 'blue', 'green'}, 1};
%   s.aSubdialog.aPopupmenu = {{10 20 30}, 1};
%   guido(s, [mfilename ' demo'])

% Note: the present m-file can be renamed without
%  requiring any changes to the Matlab code itself.
 
% svn $Id$
%=======================================================================
% Copyright (C) 1999 Dr. Charles R. Denham, ZYDECO.
%  All Rights Reserved.
%   Disclosure without explicit written consent from the
%    copyright owner does not constitute publication.
%=======================================================================
 
% Version of 15-Dec-1999 13:39:54.
% Updated    25-Jan-2000 11:47:42.

persistent CURRENT_STRUCT
persistent CURRENT_FIGURES

LF = char(10);
CR = char(13);
CRLF = [CR LF];
BABY_BLUE = [8 10 10]/10;

if nargout > 0, theResult = []; end

% Launch a demo if this is not a callback.

if nargin < 1 & isempty(gcbo)
	help(mfilename)
	ans = feval(mfilename, 'demo');
	if nargout > 0
		theResult = CURRENT_STRUCT;
	else
		assignin('caller', 'ans', CURRENT_STRUCT)
	end
	return
end

% Do the demo.

if nargin > 0 & isequal(theStruct, 'demo')
	if nargin < 2, theFigureName = [mfilename ' demo']; end
	s = [];
	s.help = help(mfilename);
	s.anEdit = 'some text';
	s.aNumber = pi;
	s.aCheckbox = {'checkbox' 0};
	s.aRadiobutton = {'radiobutton' 0};
	s.aPopupmenu = {{'red', 'blue', 'green'}, 1};
	s.aSubdialog.aPopupmenu = {{10 20 30}, 1};
%	s.aSubdialog.bSubdialog.bPopupmenu = {{100 200 300}, 1};
	disp(s)

	assign([mfilename '_demo'], s)
	result = feval(mfilename, eval([mfilename '_demo']), '', 0);
	assign([mfilename '_result'], result)
	
	if ~isempty(eval([mfilename '_result']))
		result = feval(mfilename, eval([mfilename '_result']), '', 0);
		assign([mfilename '_result'], result)
	end
	if nargout > 0
		assign('theResult', [mfilename '_result'])
	end
	return
end

if nargin > 0 & ischar(theStruct)
	theCommand = lower(theStruct);
end

if 0 & nargin > 0 & ischar(theStruct)
	switch theCommand
	case '_keypress_'
		theKey = get(gcbf, 'CurrentCharacter');
		if any(abs(theKey) < 32)
			theCommand = '_okay_';
		end
	end
end

% Dismiss dialogs in sequence.

if nargin > 0 & ischar(theStruct)
	switch theCommand
	case {'_okay_', '_cancel_'}
		CURRENT_FIGURES(~ishandle(CURRENT_FIGURES)) = [];
		if ~isempty(CURRENT_FIGURES) & gcbf ~= CURRENT_FIGURES(end)
			figure(CURRENT_FIGURES(end))
			h = helpdlg('Dismiss dialogs in sequence.', 'Please...');
			tic
			while ishandle(h) & toc < 3
				drawnow
			end
			if ishandle(h), delete(h), end
			return
		end
	end
end

% Process a command: okay or cancel.

if nargin > 0 & ischar(theStruct)
	switch theCommand
	case '_okay_'
		CURRENT_STRUCT = get(gcbf, 'UserData');
		theControls = findobj(gcbf, 'Type', 'uicontrol', 'Tag', mfilename);
		theControls = sort(theControls);
		for i = 2:2:length(theControls)
			theField = get(theControls(i-1), 'String');
			theControl = theControls(i);
			theStyle = get(theControl, 'Style');
			switch theStyle
			case 'pushbutton'
				okay = 0;
			case 'checkbox'
				theValue = get(theControl, 'Value');
				okay = 1;
			case 'radiobutton'
				theValue = get(theControl, 'Value');
				okay = 1;
			case 'edit'
				if isequal(theField, 'help')
					okay = 0;
				else
					theValue = eval(get(theControl, 'String'));
					okay = 1;
				end
			case 'popupmenu'
				theValue = get(theControl, 'Value');
				okay = 1;
			otherwise
				okay = 0;
				warning([' ## Unsupported control style: ' theStyle])
			end
			if okay
				CURRENT_STRUCT = setinfo(CURRENT_STRUCT, theField, theValue);
			end
		end
		if ~isempty(CURRENT_FIGURES) & any(CURRENT_FIGURES == gcbf)
			CURRENT_FIGURES(CURRENT_FIGURES == gcbf) = [];
		end
		delete(gcbf)
	case '_cancel_'
		CURRENT_STRUCT = [];
		if ~isempty(CURRENT_FIGURES) & any(CURRENT_FIGURES == gcbf)
			CURRENT_FIGURES(CURRENT_FIGURES == gcbf) = [];
		end
		delete(gcbf)
	otherwise
		disp([' ## Unknown command: ' theCommand])
	end
	return
end

% Build the dialog if it does not already exist.

if nargin > 0 & isstruct(theStruct)
	theStruct = setinfo(theStruct);
	if nargin < 2 | isempty(theFigureName)
		theFigureName = inputname(1);
		if isempty(theFigureName)
			theFigureName = 'unnamed';
		end
	end
	if isempty(theFigureName), theFigureName = mfilename; end
	
	f = findobj('Type', 'figure', 'Name', theFigureName, 'Tag', mfilename);
	if any(f)
		figure(f)
		return
	end
	
	if nargin < 3, isModal = 0; end
	if nargin < 4, thePrecision = 4; end
	
	if all(isModal(:)) & ~any(findstr(computer, 'PCWIN'))
		theWindowStyle = 'modal';
	else
		theWindowStyle = 'normal';
	end
	theFigures = findobj('Type', 'figure', 'Tag', mfilename);
	theFigure = figure( ...
						'Name', theFigureName, ...
						'WindowStyle', theWindowStyle, ...
						'Visible', 'off', ...
						'KeyPressFcn', [mfilename ' _keypress_'], ...
						'CloseRequestFcn', [mfilename ' _cancel_'], ...
						'Tag', mfilename);
	if any(theFigures)
		pos = get(theFigures(1), 'Position');
		left = pos(1);
		top = pos(2) + pos(4);
		for i = 2:length(theFigures)
			p = get(theFigures(i), 'Position');
			left = max(pos(1), p(1));
			top = min(top, p(2) + p(4));
		end
		thePosition = get(theFigure, 'Position');
		thePosition(1) = left + 20;
		thePosition(2) = top - thePosition(4) - 20;
		set(theFigure, 'Position', thePosition)
	end
	theFrame = uicontrol( ...
						'Style', 'frame', ...
						'Units', 'normalized', ...
						'Position', [0 0 1 1], ...
						'BackgroundColor', BABY_BLUE);
	theControls = [];
	theStruct = setinfo(theStruct);   % Canonical form.
	theFields = fieldnames(theStruct);
	for i = 1:length(theFields)
		theField = theFields{i};
		theValue = getfield(theStruct, theField);
		switch class(theValue)
		case 'cell'
			switch class(theValue{1})
			case 'cell'
				if length(theValue) > 1
					theSetting = theValue{2};
				else
					theSetting = 1;
					theValue = {theValue, theSetting};
					theStruct = setfield(theStruct, theField, theValue);
				end
				theStyle = 'popupmenu';
			case 'char'
				switch theValue{1}
				case 'checkbox'
					theStyle = 'checkbox';
					if length(theValue) > 1
						theSetting = theValue{2};
					else
						theSetting = 0;
					end
				case 'radiobutton'
					theStyle = 'radiobutton';
					if length(theValue) > 1
						theSetting = theValue{2};
					else
						theSetting = 0;
					end
				otherwise
					error([' ## Incompatible control style: ' theValue{1}])
				end
			end
			theControls(end+1) = uicontrol( ...
						'Style', 'text', ...
						'String', theField);
			theControls(end+1) = uicontrol( ...
						'Style', theStyle, ...
						'String', theValue{1}, ...
						'Value', theSetting);
		case 'char'
			f = findstr(theField, 'help_');
			if ~any(f) | f(1) ~= 1
				theControls(end+1) = uicontrol( ...
							'Style', 'text', ...
							'String', theField);
				theControls(end+1) = uicontrol( ...
							'Style', 'edit', ...
							'Max', 1000, ...
							'String', ['''' theValue '''']);
			else
				theHintName = [theField ' ' theFigureName];
				theCallback = ...
					['hint(get(gcbo, ''UserData''), ''' theHintName ''')'];
				theControls(end+1) = uicontrol( ...
							'Style', 'text', ...
							'String', theField);
				theControls(end+1) = uicontrol( ...
							'Style', 'pushbutton', ...
							'Callback', theCallback, ...
							'UserData', theValue, ...
							'String', 'Help...');
			end
		case 'double'
			theControls(end+1) = uicontrol( ...
						'Style', 'text', ...
						'String', theField);
			theControls(end+1) = uicontrol( ...
						'Style', 'edit', ...
						'Max', 1000, ...
						'String', mat2str(theValue, thePrecision));
		case 'struct'
			theCallback = ...
				[mfilename '(get(gcbo, ''UserData''), ''' theField ''')'];
			theControls(end+1) = uicontrol( ...
						'Style', 'text', ...
						'String', theField);
			theControls(end+1) = uicontrol( ...
						'Style', 'pushbutton', ...
						'String', 'More...', ...
						'Callback', theCallback, ...
						'UserData', theValue);
		otherwise
			disp(class(theValue))
			error([' ## Incompatible data type. ' class(theValue)])
		end
	end
	set(theControls(1:2:end), ...
			'HorizontalAlignment', 'right', ...
			'BackgroundColor', BABY_BLUE)
	theControls(end+1) = uicontrol( ...
			'Style', 'pushbutton', ...
			'String', 'Cancel', ...
			'BackgroundColor', [10 2 2]/10, ...
			'Callback', [mfilename ' _cancel_']);
	theControls(end+1) = uicontrol( ...
			'Style', 'pushbutton', ...
			'String', 'Okay', ...
			'BackgroundColor', [2 10 2]/10, ...
			'Callback', [mfilename ' _okay_']);
	set(theControls, 'Tag', mfilename)
	theLayout = [];
	for i = 1:length(theControls)/2
		theLayout = [(1:2:length(theControls));
					(2:2:length(theControls))].';
	end
	theLayout = theLayout(:, [1 2 2]);
	uilayout(theControls, theLayout, [1 1 18 18]/20)
	pos = get(0, 'DefaultUIControlPosition');
	width = pos(3) * 6;
	height = 0.5 * length(theControls) * pos(4) * 20 / 15;
	thePosition = get(theFigure, 'Position');
	thePosition(2) = thePosition(2) + thePosition(4) - height;
	thePosition(3) = width;
	thePosition(4) = height;
	set(theFigure, 'Position', thePosition, 'Visible', 'on', ...
					'UserData', theStruct)
	if any(CURRENT_FIGURES)
		CURRENT_FIGURES(~ishandle(CURRENT_FIGURES)) = [];
	end
	CURRENT_FIGURES(end+1) = theFigure;
	
% Wait here for the new figure to be deleted.  By then, it will
%  already have placed its contents in the persistent CURRENT_STRUCT
%  item.

	waitfor(theFigure)
	
% Now get the new info from CURRENT_STRUCT and update.
	
	theNewStruct = CURRENT_STRUCT;
	
	theFieldName = theFigureName;
	
	if length(theFigures) < 1
		theStruct = theNewStruct;
	else
		theFigure = gcf;
		theStruct = get(theFigure, 'UserData');
		if ~isfield(theStruct, theFieldName) & ~isequal(theField, 'help')
			disp([' ## No such field: ' theField])
		elseif ~isempty(theNewStruct)
			theStruct = setfield(theStruct, theFieldName, theNewStruct);
			set(theFigure, 'UserData', theStruct)
			f = findobj(theFigure, 'Type', 'uicontrol', 'Tag', mfilename);
			f = sort(f);
			for i = 2:2:length(f)
				theField = get(f(i-1), 'String');
				if isequal(theField, theFieldName)
					set(f(i), 'UserData', theNewStruct)
					break
				end
			end
		end
	end
	
	if nargout > 0
		theResult = theStruct;
	else
		assignin('caller', 'ans', theStruct);
	end
	
end

% ---------- assign ----------%

function assign(theName, theValue)

% assign -- Assign a value to a name.
%  assign('theName', theValue) assigns theValue
%   to 'theName' in the caller's workspace. It
%   avoids the need to construct an explicit
%   assignment statement to be eval-ed.
 
% Copyright (C) 1997 Dr. Charles R. Denham, ZYDECO.
%  All Rights Reserved.
%   Disclosure without explicit written consent from the
%    copyright owner does not constitute publication.
 
% Version of 28-May-1998 00:43:58.

if nargin < 2, help(mfilename), return, end

% The following scheme permits the assignment
%  of items that have complicated subscripts,
%  such as "a{1}(2).b{3}.c = pi".

hasAns = (evalin('caller', 'exist(''ans'', ''var'')') == 1);
if hasAns
	ans = evalin('caller', 'ans');   % Save.
end
assignin('caller', 'ans', theValue)
evalin('caller', [theName ' = ans;'])
evalin('caller', 'clear(''ans'')')
if hasAns
	assignin('caller', 'ans', ans)   % Restore.
end

% ---------- getinfo ----------%

function [theResult, isOkay] = getinfo(theInfo, theField)

% getinfo -- Get field value from an "Info" struct.
%  getinfo(theInfo, 'theField') returns the current
%   value of 'theField' in theInfo, a struct that
%   is compatible with the "uigetinfo" function.
%   Non-existent fields return the empty-matrix.
%  [theResult, isOkay] = ... returns isOkay = 0
%   if an error occurred; otherwise, non-zero.
%  getinfo(theInfo) returns a struct containing
%   the fields and current selections of theInfo.
 
% Copyright (C) 1997 Dr. Charles R. Denham, ZYDECO.
%  All Rights Reserved.
%   Disclosure without explicit written consent from the
%    copyright owner does not constitute publication.
 
% Version of 07-Feb-1998 09:45:56.
% Updated    15-Dec-1999 23:49:02.

if nargout > 0, theResult = []; end
if nargin < 1, help(mfilename), return, end

% Get current selections, recursively if needed.

if nargin < 2
	result = theInfo;
	theFields = fieldnames(result);
	isOkay = 1;
	for i = 1:length(theFields)
		[theValue, okay] = getinfo(theInfo, theFields{i});
		isOkay = isOkay & okay;
		switch class(theValue)
		case 'struct'
			[res, okay] = getinfo(theValue);
			isOkay = isOkay & okay;
			result = setfield(result, theFields{i}, res);
		otherwise
			result = setfield(result, theFields{i}, theValue);
		end
	end
	if nargout > 0
		theResult = result;
	else
		disp(result)
	end
	return
end

theValue = [];

isOkay = 1;
eval('theValue = getfield(theInfo, theField);', 'isOkay = 0;');

result = theValue;

if all(isOkay)
    switch class(theValue)
    case 'cell'
        if isequal(theValue{1}, 'checkbox') || ...
				isequal(theValue{1}, 'radiobutton')
			if length(theValue) < 2, theValue{2} = 0; end
            result = theValue{2};
        else
			if ~iscell(theValue{1}), theValue = {theValue{1}}; end
			if length(theValue) < 2, theValue{2} = 1; end
            result = theValue{1}{theValue{2}};
        end
    otherwise
        result = theValue;
    end
end

if nargout > 0
    theResult = result;
else
    disp(result)
end

% ---------- setinfo ----------%

function [theResult, isOkay] = setinfo(theInfo, theField, theValue)

% setinfo -- Set field value in an "Info" struct.
%  setinfo(theInfo, 'theField', theValue) updates
%   'theField' to theValue in theInfo, a struct
%   that is compatible with the "uigetinfo" function.
%   If 'theField' does not exist, it will be created
%   to receive theValue.
%  [theResult, isOkay] = ... returns isOkay = 0
%   if an error occurred; otherwise, non-zero.
%  setinfo(theInfo, 'theField') invokes "getinfo".
%  setinfo(theInfo) adjusts theInfo to canonical
%   form.
 
% Copyright (C) 1997 Dr. Charles R. Denham, ZYDECO.
%  All Rights Reserved.
%   Disclosure without explicit written consent from the
%    copyright owner does not constitute publication.
 
% Version of 07-Feb-1998 09:45:56.
% Updated    16-Dec-1999 01:20:36.

if nargout > 0, theResult = []; end
if nargin < 1, help(mfilename), return, end

% Put all fields into canonical form.

if nargin < 2
	result = theInfo;
	theFields = fieldnames(result);
	isOkay = 1;
	for i = 1:length(theFields)
		theValue = getfield(result, theFields{i});
		switch class(theValue)
		case 'struct'
			[res, okay] = setinfo(theValue);
			isOkay = isOkay & okay;
			result = setfield(result, theFields{i}, res);
		case 'cell'
			if length(theValue) == 1 && ...
				(isequal(theValue{1}, 'radiobutton') || ...
				isequal(theValue{1}, 'checkbox'))
				theValue{2} = 0;
				result = setfield(result, theFields{i}, theValue);
			elseif length(theValue) == 1 && iscell(theValue{1})
				theValue{2} = 1;
				result = setfield(result, theFields{i}, theValue);
			elseif ~isequal(theValue{1}, 'radiobutton') && ...
					~isequal(theValue{1}, 'checkbox')
				theValue{1} = {theValue{1}};
				theValue{2} = 1;
			else
				theValue;
			end
		otherwise
		end
	end
	if nargout > 0, theResult = result; end
	return
end

[theVal, isOkay] = getinfo(theInfo, theField);

if nargin == 2
    if nargout > 0
        theResult = theVal;
    else
        disp(theVal)
    end
    return
end

result = theInfo;

if ~all(isOkay)   % Create a new field.
    isOkay = 1;
    eval('result = setfield(theInfo, theField, theValue);', 'isOkay = 0;');
else   % Update an existing field.
	isokay = 1;
	eval('theVal = getfield(theInfo, theField);', 'isOkay = 0;');
	if ~isOkay, theVal = []; end
    switch class(theVal)
	case 'cell'
        if isequal(theVal{1}, 'checkbox') | isequal(theVal{1}, 'radiobutton')
			% Do nothing.
		elseif ~iscell(theVal{1})
			theVal{1} = {theVal{1}};
			theVal{2} = 1;
		end
	end
    switch class(theVal)
    case 'cell'
        if isequal(theVal{1}, 'checkbox') | isequal(theVal{1}, 'radiobutton')
            theVal{2} = any(any(theValue));
        else
			switch class(theValue)
			case 'double'
				if theValue > 0 & theValue <= length(theVal{1})
					theVal{2} = theValue;
				end
			otherwise
            	flag = 0;
	            for i = 1:length(theVal{1})
	                if isequal(theVal{1}{i}, theValue)
	                    theVal{2} = i;
	                    flag = 1
	                end
	            end
	            if ~any(flag)   % Append.
	                theVal{1} = [theVal(:); {theValue}];
	                theVal{2} = length(theVal{1});
	            end
			end
% else
% theVal{1} = [{theValue}; theVal(:)];
        end
    otherwise
        theVal = theValue;
    end
    isOkay = 1;
    eval('result = setfield(theInfo, theField, theVal);', 'isOkay = 0;');
end

if nargout > 0
    theResult = result;
else
    disp(result)
end

% ---------- uilayout ----------%

function theResult = uilayout(theControls, theLayout, thePosition)

% uilayout -- Layout for ui controls.
%  uilayout(theControls, theLayout) positions theControls
%   according to theLayout, an array whose entries, taken
%   in sorted order, define the rectangular extents occupied
%   by each control.  TheLayout defaults to a simple vertical
%   arrangement of theControls.  A one-percent margin is
%   imposed between controls.  To define a layout region
%   containing no control, use Inf.
%  uilayout(..., thePosition) confines the controls to the
%   given normalized position of the figure.  This syntax
%   is useful for embedding controls within a frame.
%  uilayout (no argument) demonstrates itself.
 
% Copyright (C) 1997 Dr. Charles R. Denham, ZYDECO.
%  All Rights Reserved.
%   Disclosure without explicit written consent from the
%    copyright owner does not constitute publication.
 
% Version of 18-Apr-1997 08:07:54.
% Updated    27-Dec-1999 06:03:57.

if nargin < 1, theControls = 'demo'; help(mfilename), end

if strcmp(theControls, 'demo')
   theLayout = [1 2;
                3 4;
                5 Inf;
                5 6;
                5 Inf;
                7 8;
                9 10;
                11 12;
                13 14];
   [m, n] = size(theLayout);
   thePos = get(0, 'DefaultUIControlPosition');
   theSize = [n+2 m+2] .* thePos(3:4);
   theFigure = figure('Name', 'UILayout', ...
                      'NumberTitle', 'off', ...
                      'Resize', 'off', ...
                      'Units', 'pixels');
   thePos = get(theFigure, 'Position');
   theTop = thePos(2) + thePos(4);
   thePos = thePos .* [1 1 0 0] + [0 0 theSize];
   thePos(2) = theTop - (thePos(2) + thePos(4));
   set(theFigure, 'Position', thePos);
   theFrame = uicontrol('Style', 'frame', ...
                        'Units', 'normalized', ...
                        'Position', [0 0 1 1], ...
                        'BackgroundColor', [0.5 1 1]);
   theStyles = {'checkbox'; 'text'; ...
                'edit'; 'text'; ...
                'listbox'; 'text'; ...
                'popupmenu'; 'text'; ...
                'pushbutton'; 'text'; ...
                'radiobutton'; 'text'; ...
                'text'; 'text'};
   theStrings = {'Anchovies?', '<-- CheckBox --', ...
                 'Hello World!', '<-- Edit --', ...
                 {'Now', 'Is', 'The' 'Time' 'For' 'All' 'Good', ...
                  'Men', 'To', 'Come' 'To' 'The' 'Aid' 'Of', ...
                  'Their' 'Country'}, ...
                 '<-- ListBox --', ...
                 {'Cheetah', 'Leopard', 'Lion', 'Tiger', 'Wildcat'}, ...
                 '<-- PopupMenu --', ...
                 'Okay', '<-- PushButton --', ...
                 'Cream?', '<-- RadioButton --', ...
                 'UILayout', '<-- Text --'};
   theControls = zeros(size(theStyles));
   for i = 1:length(theStyles)
      theControls(i) = uicontrol('Style', theStyles{i}, ...
                                 'String', theStrings{i}, ...
                                 'Callback', ...
                                 'disp(int2str(get(gcbo, ''Value'')))');
   end
   set(theControls(1:2:length(theControls)), 'BackGroundColor', [1 1 0.5])
   set(theControls(2:2:length(theControls)), 'BackGroundColor', [0.5 1 1])
   thePosition = [1 1 98 98] ./ 100;
   uilayout(theControls, theLayout, thePosition)
   set(theFrame, 'UserData', theControls)
   theStyles, theLayout, thePosition
   if nargout > 0, theResult = theFrame; end
   return
end

if nargin < 2, theLayout = (1:length(theControls)).'; end
if nargin < 3, thePosition = [0 0 1 1]; end

a = theLayout(:);
a = a(isfinite(a));
a = sort(a);
a(diff(a) == 0) = [];

b = zeros(size(theLayout));

for k = 1:length(a)
   b(theLayout == a(k)) = k;
end

[m, n] = size(theLayout);

set(theControls, 'Units', 'Normalized')
theMargin = [1 1 -2 -2] ./ 100;
for k = 1:min(length(theControls), length(a))
   [i, j] = find(b == k);
   xmin = (min(j) - 1) ./ n;
   xmax = max(j) ./ n;
   ymin = 1 - max(i) ./ m;
   ymax = 1 - (min(i) - 1) ./ m;
   thePos = [xmin ymin (xmax-xmin) (ymax-ymin)] + theMargin;
   thePos = thePos .* thePosition([3 4 3 4]);
   thePos(1:2) = thePos(1:2) + thePosition(1:2);
   set(theControls(k), 'Position', thePos);
end

if nargout > 0, theResult = theControls; end

% ---------- hint ----------%

function theResult = hint(theText, theFigureName)

% hint -- Post a message.
%  hint('theText', 'theFigureName') posts 'theText' in the
%   figure named 'theFigureName', which will be created if
%   it does not already exist.  The text can be a string
%   separated by newlines, or a cell of strings.
%  hint('demo') demonstrates itself.
 
% Copyright (C) 1999 Dr. Charles R. Denham, ZYDECO.
%  All Rights Reserved.
%   Disclosure without explicit written consent from the
%    copyright owner does not constitute publication.
 
% Version of 17-Dec-1999 11:31:03.
% Updated    04-Jan-2000 13:40:41.

if nargin < 1, theText = 'demo'; end

if isequal(theText, 'demo')
	hint(help(mfilename), 'Hint Hint')
	return
end

CR = char(13);
LF = char(10);
CRLF = [CR LF];

if nargout > 0, theResult = []; end
if nargin < 1
	help(mfilename)
	theText = help(mfilename);
	for i = 1:3
		theText = [theText theText];
	end
end
if nargin < 2, theFigureName = 'Help'; end

if min(size(theText)) == 1
	theText = theText(:).';
end

if ischar(theText) & size(theText, 1) == 1
	theText = strrep(theText, CRLF, CR);
	theText = strrep(theText, LF, CR);
	if theText(end) ~= CR, theText(end+1) = CR; end
	if theText(1) ~= CR, theText = [CR theText]; end
	f = find(theText == CR);
	if any(f)
		s = cell(length(f)-1, 1);
		for i = 1:length(f)-1
			if f(i)+1 < f(i+1)-1
				s{i} = theText(f(i)+1:f(i+1)-1);
			else
				s{i} = ' ';
			end
		end
		theText = s;
	end
end

theLineCount = size(theText, 1);

theFigure = findobj( ...
					'Type', 'figure', ...
					'Name', theFigureName, ...
					'Tag', [mfilename ' hint'] ...
					);
if isempty(theFigure)
	theFigures = findobj('Type', 'figure');
	if any(theFigures)
		thePosition = get(gcf, 'Position');
		left = thePosition(1);
		top = thePosition(2) + thePosition(4);
	end
	theFigure = figure( ...
					'Name', theFigureName, ...
					'Tag', mfilename, ...
					'Visible', 'off' ...
					);
	if any(theFigures)
		set(theFigure, 'Position', thePosition + [20 -20 0 0])
	end
end

delete(get(theFigure, 'Children'))

theColor = [9 9 9]/10;
theFontName = 'Courier';
if any(findstr(lower(computer), 'mac'))
	theFontName = 'Monaco';
end
theFontSize = 12;
theListBox = uicontrol( ...
				theFigure, ...
				'Style', 'listbox', ...
				'String', theText, ...
				'FontName', theFontName, ...
				'FontSize', theFontSize, ...
				'HorizontalAlignment', 'left', ...
				'BackgroundColor', theColor ...
				);

uilayout(theListBox, 1)

set(theFigure, 'Visible', 'on')
figure(theFigure)

if nargout > 0
	theResult = theFigure;
end
