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
 
% svn $Id$
%=======================================================================
% Copyright (C) 1997 Dr. Charles R. Denham, ZYDECO.
%  All Rights Reserved.
%   Disclosure without explicit written consent from the
%    copyright owner does not constitute publication.
%=======================================================================
 
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
			if length(theValue) == 1 & ...
				(isequal(theValue{1}, 'radiobutton') | ...
				isequal(theValue{1}, 'checkbox'))
				theValue{2} = 0;
				result = setfield(result, theFields{i}, theValue);
			elseif length(theValue) == 1 & iscell(theValue{1})
				theValue{2} = 1;
				result = setfield(result, theFields{i}, theValue);
			elseif ~isequal(theValue{1}, 'radiobutton') & ...
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
