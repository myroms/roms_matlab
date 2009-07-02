function varargout = inherit(method, self, varargin)

% inherit -- Inherit a superclass method.
%  varargout = inherit('method', self, varargin) calls
%   the superclass 'method' of self, an object, with the
%   given input and output arguments.  The routine
%   climbs the inheritance tree if needed.

% svn $Id$
%=======================================================================
% Copyright (C) 1996 Dr. Charles R. Denham, ZYDECO.
%  All Rights Reserved.
%   Disclosure without written consent from the
%    copyright owner does not constitute publication.
%=======================================================================

if nargin < 2, help inherit, return, end

if ~isobject(self), return, end

theSuperObject = super(self);
varargin = [cell(1, 0); varargin];

while isobject(theSuperObject)
   varargin{1} = theSuperObject;
   failure = 0;
   eval(vargstr(method, varargin, varargout), 'failure = 1;')
   if ~failure, break; end
   theSuperObject = super(theSuperObject);
end
