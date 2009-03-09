function theResult = verbose(self, isVerbose)

% verbose -- Set/get verbosity flag.
%  verbose(self, isVerbose) sets the verbosity
%   flag to isVerbose.
%  verbose(self) returns the verbosity flag.
 
% svn $Id$
%=======================================================================
% Copyright (C) 2000 Dr. Charles R. Denham, ZYDECO.
%  All Rights Reserved.
%   Disclosure without explicit written consent from the
%    copyright owner does not constitute publication.
%=======================================================================
 
% Version of 05-Feb-2000 13:02:56.
% Updated    05-Feb-2000 13:02:56.

if nargin > 1
	if ischar(isVerbose)
		isVerbose = isequal(isVerbose, 'on');
	end
	isVerbose = ~~isVerbose;
	psset(self, 'itIsVerbose', isVerbose);
end

isVerbose = psget(self, 'itIsVerbose');

if nargout > 0
	theResult = isVerbose;
elseif isVerbose
	disp(' ## Verbose on')
else
	disp(' ## Verbose off')
end
