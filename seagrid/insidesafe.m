function theResult = insidesafe(x, y, px, py, theChunkSize, theVerboseFlag)

% insidesafe -- Points in polygon via "inpolygon".
%  insidesafe('demo') demonstrates itself with 100 random points.
%  insidesafe(N) demonstrates itself with N random points (default = 100).
%  insidesafe(x, y, px, py) calls "inpolygon" to determine
%   which (x, y) points are inside the polygon with vertices
%   (px, py).  Result is 1 for inside, 0 for outside, and 0.5
%   for on-the-line.
%  insidesafe(..., theChunkSize) starts in steps of theChunkSize,
%   which is then reduced whenever an "out-of-memory" error is
%   encountered.  The default value is 25.
%  insidesafe(..., theChunkSize, theVerboseFlag) displays progress
%   information if theVerboseFlag is TRUE.
 
% svn $Id$
%=======================================================================
% Copyright (C) 1999 Dr. Charles R. Denham, ZYDECO.
%  All Rights Reserved.
%   Disclosure without explicit written consent from the
%    copyright owner does not constitute publication.
%=======================================================================
 
% Version of 28-Apr-1999 20:06:47.
% Updated    16-Mar-2001 11:42:11.

% N.B. -- Another scheme would be to triangulate the polygon,
%  then to use "tsearch" to find points inside the convex hull,
%  which would need to be trimmed to eliminate concave regions
%  at the convex-hull boundary.  Not a simple game.  The first
%  part is very fast, but the culling is dreadfully slow.

DEFAULT_CHUNK_SIZE = 25;  % <== Adjust as needed.
DEFAULT_VERBOSE_FLAG = ~~0;   % Verbosity.

% Demonstration.

if nargin < 1, x = 'demo'; help(mfilename), end
if isequal(x, 'demo'), x = 100; end
if ischar(x), x = eval(x); end
if length(x) == 1
	px = [0 1 1 0 0]/2;
	py = [0 0 1 1 0]/2;
	n = x;
	x = rand(1, n);
	y = rand(size(x));
	tic
	result = insidesafe(x, y, px, py, DEFAULT_CHUNK_SIZE, 0);
	t = toc;
	disp([' ## Elapsed time: ' num2str(t) ' s.'])
	k = (result == 1);
	h = plot(px, py, x(k), y(k), 'go', x(~k), y(~k), 'r+');
	legend(h, 'polygon', 'inside', 'outside')
	figure(gcf)
	set(gcf, 'Name', [mfilename ' ' int2str(n)])
	return
end

if nargin < 4, help(mfilename), return, end
if nargin < 5, theChunkSize = DEFAULT_CHUNK_SIZE; end
if nargin < 6, theVerboseFlag = DEFAULT_VERBOSE_FLAG; end

px = px(:); py = py(:);
if length(px) ~= length(py)
	error(' ## Polygon px and py must be same length.')
end

if length(x) == 1, x = x + zeros(size(y)); end
if length(y) == 1, y = y + zeros(size(x)); end

result = logical(zeros(size(x)));
kmax = prod(size(result));

okay = 1;

tic

milestone = kmax - rem(kmax, theChunkSize);

k = 0;
while k < kmax
	if theVerboseFlag
		remaining = kmax - k;
		if k > 0 & remaining < milestone
			time = remaining * toc / k;
			disp([' ## Remaining: ' int2str(remaining) ' points; ' int2str(time) ' seconds.'])
			milestone = milestone - 4*theChunkSize;
		end
	end
	i = k+1:min(k+theChunkSize, kmax);
	try
		result(i) = inpolygon(x(i), y(i), px, py);
		k = k + theChunkSize;
	catch
		if findstr(lower(lasterr), 'memory')
			lasterr('')
			theChunkSize = floor(theChunkSize / 2);
			disp([' ## Chunk-size reduced to: ' int2str(theChunkSize)])
			if theChunkSize < 1
				error(' ## Chunk-size reduced to zero.')
			end
		else
			okay = 0;
			break
		end
	end
end

if ~okay
	if isempty(lasterr)
		error([' ## Interrupted.'])
	else
		error([' ## ' lasterr])
	end
end

if nargout > 0
	theResult = result;
else
	disp(result)
end
