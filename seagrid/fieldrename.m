function theResult = fieldrename(theStruct, varargin)

% fieldrename -- Rename fields in a "struct".
%  fieldrename(theStruct, 'theFieldName', 'theNewName') renames
%   'theFieldName' to 'theNewName' in theStruct.  The arguments
%   can be entered on the command line without parentheses, in
%   which case, 'theStruct' will be treated as the name of the
%   actual struct in the caller's workspace.  More than one
%   fieldname can be changed during a single call by appending
%   additional fieldname/newname pairs to the argument list.
%   An empty new-name ('') causes that particular field to be
%   eliminated from the structure.
 
% svn $Id$
%=======================================================================
% Copyright (C) 1999 Dr. Charles R. Denham, ZYDECO.
%  All Rights Reserved.
%   Disclosure without explicit written consent from the
%    copyright owner does not constitute publication.
%=======================================================================
 
% Version of 13-Sep-1999 11:45:49.
% Updated    20-Dec-1999 15:56:48.

if nargin < 3, return, end

theStructName = inputname(1);

if ischar(theStruct)
	theStructName = theStruct;
	theStruct = evalin('caller', theStructName);
end

update = 0;

f = fieldnames(theStruct);
g = f;
for j = 1:2:length(varargin)
	theFieldName = varargin{j};
	theNewName = varargin{j+1};
	if isequal(theNewName, '''''') | isequal(theNewName, '[]')
		theNewName = '';
	end
	for i = length(g):-1:1
		if isequal(g{i}, theFieldName)
			if isempty(theNewName)
				g(i) = [];
			else
				g{i} = theNewName;
			end
			update = 1;
			break
		end
	end
end

if update & length(g) > 0
	result = cell(1, 2*length(g));
	k = 0;
	for i = 1:length(g)
		result{k+1} = g{i};
		result{k+2} = cell(size(theStruct));
		for j = 1:prod(size(theStruct))
			theValue = getfield(theStruct(j), f{i});
			result{k+2}{j} = theValue;
		end
		k = k+2;
	end
	result = struct(result{:});
else
	result = [];
end

if nargout > 0
	theResult = result;
else
	assignin('caller', theStructName, result)
end
