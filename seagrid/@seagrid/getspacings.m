function theResult = getspacings(self)

% seagrid/getspacings -- Get default spacings.
%  getspacings(self) invokes a dialog to get functions
%   of 's' that describe the default positions of grid
%   lines along edges #1 and #2, on behalf of self, a
%   "seagrid" object.  The functions are "evaled", using
%   vector s = linspace(0, 1, N), where N is the number
%   of points along the respective edge.  Thus, 's' alone
%   will provide a uniform spacing; 's.^2' will produce
%   quadratic spacing, i.e., a linearly-decreasing
%   density of grid lines; etc.  Note that the function
%   should use scalar operators.  The result will be
%   sorted and normalized to fit the range [0:1].
%
%   Since the result is sorted, the spacing functions
%   need not be monotonically increasing over its range.
%
%   If the entry is a vector surrounded by [...], it will
%   be treated as a set of "way-points", to be interpolated
%   sensibly to fit the actual number of existing grid lines.
%   (NOT YET IMPLEMENTED)
 
% svn $Id$
%=======================================================================
% Copyright (C) 1999 Dr. Charles R. Denham, ZYDECO.
%  All Rights Reserved.
%   Disclosure without explicit written consent from the
%    copyright owner does not constitute publication.
%=======================================================================
 
% Version of 03-Aug-1999 10:37:40.
% Updated    07-Mar-2000 22:42:48.

if nargout > 0, theResult = self; end
if nargin < 1, help(mfilename), return, end

theGridSize = psget(self, 'itsGridSize');

p = psget(self, 'itsSpacings');
q = psget(self, 'itsDefaultSpacings');
e = psget(self, 'itsSpacedEdges');

if isempty(p)
	p = cell(1, length(theGridSize));
	for k = 1:length(theGridSize)
		p{k} = linspace(0, 1, theGridSize(k)+1);
	end
	psset(self, 'itsSpacings', p)
end

if isempty(q)
	q = {'s', 's', 0};
	psset(self, 'itsDefaultSpacings', q)
end

if isempty(e), e = [1 2]; end
selector = (e == [1 2]);
selector(~selector) = 2;

d.Spacings_Edge_A = q{1};   % uigetinfo dialog structure.
d.Select_Edge_A = {{1, 3}, selector(1)};
d.Spacings_Edge_B = q{2};
d.Select_Edge_B = {{2, 4}, selector(2)};
% d.Density_Flag = {'checkbox', q{3}};
SeaGrid_Spacing_Setup = d;

if (0)
%	reply = uigetinfo(SeaGrid_Spacing_Setup);
else
	reply = guido(SeaGrid_Spacing_Setup);
end

% If the Density_Flag is checked, then we will have
%  to integrate and normalize the given expressions.
%  The flexibility is needed because some direct
%  distributions might not be easilty expressable.
%  *** Not yet done. ***

if ~isempty(reply)
	q{1} = getinfo(reply, 'Spacings_Edge_A');
	q{2} = getinfo(reply, 'Spacings_Edge_B');
	q{3} = getinfo(reply, 'Density_Flag');
	if isempty(q{3}), q{3} = 0; end
	psset(self, 'itsDefaultSpacings', q);
	e(1) = getinfo(reply, 'Select_Edge_A');
	e(2) = getinfo(reply, 'Select_Edge_B');
	psset(self, 'itsSpacedEdges', e)
	for k = 1:2
		s = linspace(0, 1, length(p{k}));
		if ~isempty(q{k})
			bad = 0;
			eval(['s = ' q{k} ';'], 'bad = 1;');
			if bad
				disp([' ## Unable to evaluate: ' q{k}])
			end
			s = sort(s);
			s = s - min(s); s = s / max(s);
		end
		p{k} = s;
	end
	psset(self, 'itsSpacings', p)
	doupdate(self, 1)
end

if nargout > 0, theResult = self; end
