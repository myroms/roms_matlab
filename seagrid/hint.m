function theResult = hint(theText, theFigureName)

% hint -- Post a message.
%  hint('theText', 'theFigureName') posts 'theText' in the
%   figure named 'theFigureName', which will be created if
%   it does not already exist.  The text can be a string
%   separated by newlines, or a cell of strings.
%  hint('demo') demonstrates itself.
 
% svn $Id$
%=======================================================================
% Copyright (C) 1999 Dr. Charles R. Denham, ZYDECO.
%  All Rights Reserved.
%   Disclosure without explicit written consent from the
%    copyright owner does not constitute publication.
%=======================================================================
 
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
					'Tag', mfilename ...
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
