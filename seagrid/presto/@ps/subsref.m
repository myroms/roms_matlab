function theResult = subsref(self, theStruct)

% ps/subsref -- Process subscripted references.
%  subsref(self, theStruct) processes subscripted
%   references on behalf of self, a "ps" object.
 
% svn $Id$
%=======================================================================
% Copyright (C) 1999 Dr. Charles R. Denham, ZYDECO.
%  All Rights Reserved.
%   Disclosure without explicit written consent from the
%    copyright owner does not constitute publication.
%=======================================================================
 
% Version of 09-Dec-1999 01:19:01.
% Updated    09-Dec-1999 01:19:01.

if nargin < 2, help(mfilename), return, end

theType = theStruct(1).type;
theSubs = theStruct(1).subs;
if ~iscell(theSubs), theSubs = {theSubs}; end
theStruct(1) = [];

switch theType
case '.'
	result = psget(self, theSubs{1});
	if length(theStruct) > 0
		result = subsref(result, theStruct);
	end
otherwise
end

theResult = result;
