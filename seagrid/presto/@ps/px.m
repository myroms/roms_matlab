function [varargout] = px(varargin)

% svn $Id$
%=======================================================================
% Copyright (C) 1999 Dr. Charles R. Denham, ZYDECO.
%  All Rights Reserved.
%   Disclosure without explicit written consent from the
%    copyright owner does not constitute publication.
%=======================================================================
 
% Version of 17-Dec-1999 16:12:35.
% Updated    17-Dec-1999 16:12:35.

disp([' ## Obsolete function: ' mfilename])

varargout = cell(1, nargout);

if nargout > 0
	[varargout{:}] = ps(varargin{:});
else
	ps(varargin{:});
end
