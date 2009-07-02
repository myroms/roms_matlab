function theSuperObject = super(theObject)

% super -- Super-object of an object.
%  super(theObject) returns the super-object
%   of theObject, or [] if none exists.  The
%   super-object is the bottom-most object
%   in the struct of theObject, if any.

% svn $Id$
%=======================================================================
% Copyright (C) 1996 Dr. Charles R. Denham, ZYDECO.
%  All Rights Reserved.
%   Disclosure without written consent from the
%    copyright owner does not constitute publication.
%=======================================================================

% Version of 04-Apr-1997 16:51:36.
% Revised    02-Nov-1998 08:26:00.

if nargin < 1, help super, return, end

if isobject(theObject)
   theStruct = struct(theObject);
  else
   theStruct = theObject;
end

s = [];

f = fieldnames(theStruct);
if ~isempty(f)
   s = getfield(theStruct, f{length(f)});
   if ~isobject(s), s = []; end
end

if nargout > 0
   theSuperObject = s;
  else
   disp(s)
end
