function theResult = ps(theHandle)

% ps/ps -- Constructor for "ps" class.
%  ps(theHandle) attaches theHandle to a new "ps"
%   object, unless the handle is already associated
%   with such an object.  In either case, the object
%   itself is returned.  The default graphics entity
%   is a new figure named "PS".
%  ps(aPSObject) returns the handle associated
%   with the given "ps" object.
% Also see: "help ps_test".
 
% svn $Id$
%=======================================================================
% Copyright (C) 1999 Dr. Charles R. Denham, ZYDECO.
%  All Rights Reserved.
%   Disclosure without explicit written consent from the
%    copyright owner does not constitute publication.
%=======================================================================
 
% Version of 09-Dec-1999 00:12:19.
% Updated    09-Dec-1999 05:26:41.

if nargout < 1 & nargin < 1
	setdef(mfilename)
	return
end

if nargout > 0, theResult = []; end

if nargin < 1
	help(mfilename)
	theHandle = figure('Name', 'PS');
end

% Convert existing handle to "ps" object.

switch class(theHandle)
case 'double'
	if ishandle(theHandle)
		u = get(theHandle, 'UserData');
		if isstruct(u) & isfield(u, 'ps_Self')
			theSelf = u.ps_Self;
			if isa(theSelf, 'ps')
				self = theSelf;
				if nargout > 0
					theResult = self;
				else
					assignin('caller', 'ans', self)
				end
				return
			end
		end
	end
end

% Convert existing "ps" object to its handle.

if isa(theHandle, 'ps')
	self = theHandle;
	theHandle = self.itsHandle;
	if nargout > 0
		theResult = theHandle;
	else
		assignin('caller', 'ans', theHandle)
	end
	return
end

% Create a new "ps" object.
		
theStruct.itsHandle = [];
self = class(theStruct, 'ps');
self.itsHandle = theHandle;

% Note names of the key elements: "ps_Self" and "ps_Data".
%  These appear also in "isps()" and "ps/data()".  Make
%  sure they are used consistently throughout.

psbind(self)

if nargout > 0
	theResult = self;
else
	assignin('caller', 'ans', self)
end
