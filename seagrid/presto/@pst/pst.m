function theResult = pst(varargin)

% pst/pst -- Constructor for "pst" class.
%  pst(...) constructs a "pst" object from the given
%   arguments.  Use this class as a template for
%   deriving a customized class from "ps".
 
% svn $Id$
%=======================================================================
% Copyright (C) 1999 Dr. Charles R. Denham, ZYDECO.
%  All Rights Reserved.
%   Disclosure without explicit written consent from the
%    copyright owner does not constitute publication.
%=======================================================================
 
% Version of 05-Nov-1999 20:11:19.
% Updated    05-Nov-1999 20:11:19.

theStruct.ignore = [];

myPS = ps(figure);
theFigure = handle(myPS);

self = class(theStruct, 'pst', myPS);

u = get(theFigure, 'UserData');   % We need a "bind" routine.
u.ps_Self = self;
set(theFigure, 'UserData', u)

if nargout  0
	theResult = self;
end
