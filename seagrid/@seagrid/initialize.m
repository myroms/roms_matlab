function theResult = initialize(self, varargin)

% seagrid/initialize -- Initialization of "seagrid".
%  initialize(self, ...) initializes self, a "seagrid" object,
%   with the name/value pairs given in the argument-list.
 
% svn $Id$
%=======================================================================
% Copyright (C) 1999 Dr. Charles R. Denham, ZYDECO.
%  All Rights Reserved.
%   Disclosure without explicit written consent from the
%    copyright owner does not constitute publication.
%=======================================================================
 
% Version of 07-Apr-1999 15:47:38.
% Updated    03-Jun-1999 17:24:28.

if nargin < 1, help(mfilename), return, end

for i = 2:2:length(varargin)
	self = psset(self, varargin{i-1}, varargin{i})
end

theProjection = psget(self, 'itsProjection');
m_proj(theProjection)

getboundary(self)

if nargout > 0, theResult = self; end
