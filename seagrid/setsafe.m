function theResult = setsafe(varargin)

% setsafe -- Fault-tolerant and non-redundant "set".
%  setsafe('demo') demonstrates itself.
%  setsafe(theHandle, ...) sets only those properties
%   that would actually change.  The argument-list
%   consists of a valid handle, followed by name/value
%   pairs, followed by another handle, etc., so long
%   as the list makes sense.  In the case of multiple
%   handles, each block is set before processing the next,
%   as if several "set" calls were executed in sequence.
%   The result is logically true if some values were set
%   successfully; otherwise, it is false.
%
% Note: this routine overcomes the propensity of Matlab
%  to generate redundant screen-updates under X-Windows.
 
% svn $Id$
%=======================================================================
% Copyright (C) 2000 Dr. Charles R. Denham, ZYDECO.
%  All Rights Reserved.
%   Disclosure without explicit written consent from the
%    copyright owner does not constitute publication.
%=======================================================================
 
% Version of 11-Feb-2000 11:56:50.
% Updated    11-Feb-2000 17:02:36.

global SETSAFE_VERBOSE

if nargout > 0, theResult = 1; end
if nargin < 1, help(mfilename), varargin = {'demo'}; end

if isempty(SETSAFE_VERBOSE), SETSAFE_VERBOSE = ~~0; end

if isequal(varargin{1}, 'demo')
	SETSAFE_VERBOSE = ~~1;
	disp([' >> ' mfilename ' demo']), disp(' ')
	c{1} = 'setsafe(gcf, ''Name'', [mfilename '' demo''])';
	c{2} = 'setsafe(gcf, ''Color'', [1 0 0], ''Color'', [0 0 1])';
	disp([' >> ' c{1}])
	disp([' >> ' c{2}])
	eval(c{1})
	eval(c{2})
	SETSAFE_VERBOSE = ~~0;
	return
end

vargs = {};
mode = 'handle';
i = 0;  % Index of arguments already processed.

while i <= length(varargin)
	if i >= length(varargin)
		mode = 'set';
	else
		v = varargin{i+1};
	end
	switch mode
	case 'handle'   % Looking for a handle.
		if ishandle(v)
			if SETSAFE_VERBOSE, disp([' ## ' mode ' ' mat2str(v)]), end
			vargs = {v};
			mode = 'property';
			i = i+1;
		else
			mode = 'ignore';
		end
	case 'property'   % Looking for a property-name.
		if ischar(v)
			if SETSAFE_VERBOSE, disp([' ## ' mode ' ''' mat2str(v) '''']), end
			vargs{end+1} = v;
			mode = 'value';
			i = i+1;
		elseif ishandle(v)
			mode = 'set';
		else
			mode = 'ignore';
		end
	case 'value'   % Looking for a property-value.
		if SETSAFE_VERBOSE
			theClass = class(v);
			switch theClass
			case 'double'
				disp([' ## ' mode ' ' mat2str(v)])
			case 'char'
				disp([' ## ' mode ' ''' mat2str(v) ''''])
			otherwise
				disp([' ## ' mode])
				disp(v)
			end
		end
		vargs{end+1} = v;
		mode = 'property';
		i = i+1;
	case 'set'
		if SETSAFE_VERBOSE, disp([' ## ' mode]), end
		if nargout > 0
			theResult = nrset(vargs{:});
		else
			nrset(vargs{:});
		end
		if i == length(varargin), break; end
		vargs = {};
		mode = 'handle';
	case 'ignore'
		i = i+1;
		mode = 'set';
	end
	if i >= length(varargin)
		mode = 'set';
	end
end

% ---------- nrset ---------- %

function theResult = nrset(varargin)

% nrset -- Non-redundant "set".
%  nrset(theHandle, ...) performs non-redundant "set",
%   which attempts to set only those properties that
%   would actually change.  All of the requested
%   properties are checked before any of them are set.
%   The result is logically true if any values were
%   set successfully; otherwise, it is false.
 
% Copyright (C) 2000 Dr. Charles R. Denham, ZYDECO.
%  All Rights Reserved.
%   Disclosure without explicit written consent from the
%    copyright owner does not constitute publication.
 
% Version of 11-Feb-2000 12:57:27.
% Updated    20-Apr-2000 15:23:10.

global SETSAFE_VERBOSE

if nargout > 0, theResult = 0; end
if nargin < 1, help(mfilename), return, end

if length(varargin) < 1, return, end

theHandle = varargin{1};

k = length(varargin);
if rem(k, 2), k = k-1; end

while k > 1
	theName = varargin{k};
	if k < length(varargin)
		theValue = varargin{k+1};
		try
			redundant = ~~1;
			for i = 1:length(theHandle)
				x = get(theHandle(1), theName);
				if ~isequal(x, theValue)
					redundant = ~~0;
					break
				end
			end
			if redundant
				if any(SETSAFE_VERBOSE)
					disp([' ## redundant property ''' theName ''''])
				end
				varargin(k:k+1) = [];
			end
		catch
			disp([' ## ' mfilename ': ' lasterr])
			return
		end
	end
	k = k-2;
end

okay = (nargin > 1 & length(varargin) > 1) | ...
			(nargin == 1 & length(varargin) == 1);

if okay
	try
		if nargout > 0
			theResult = set(varargin{:});
		else
			set(varargin{:});
		end
		theResult = 1;
	catch
		disp([' ## ' mfilename ': ' lasterr])
	end
end
