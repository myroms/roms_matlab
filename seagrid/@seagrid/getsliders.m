function theResult = getsliders(self)

% seagrid/getsliders -- Get slider count.
%  getsliders(self) invokes a dialog to get the
%   number of sliders along edges #1 and #2.
 
% svn $Id$
%=======================================================================
% Copyright (C) 1999 Dr. Charles R. Denham, ZYDECO.
%  All Rights Reserved.
%   Disclosure without explicit written consent from the
%    copyright owner does not constitute publication.
%=======================================================================
 
% Version of 03-Aug-1999 10:37:40.
% Updated    10-Jan-2000 16:47:03.

if nargin < 1, help(mfilename), return, end

theGridSize = psget(self, 'itsGridSize');

p = psget(self, 'itsSpacings');
q = psget(self, 'itsDefaultSpacings');

if isempty(p)
	p = cell(1, length(theGridSize));
	for k = 1:length(theGridSize)
		p{k} = linspace(0, 1, theGridSize(k));
	end
	psset(self, 'itsSpacings', p)
end

if isempty(q)
	q = {'s', 's', 0};
	psset(self, 'itsDefaultSpacings', q)
end

d.Spacings_Edge_1 = q{1};   % uigetinfo dialog structure.
d.Spacings_Edge_2 = q{2};
d.Density_Flag = {'checkbox', q{3}};
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
	q{1} = getinfo(reply, 'Spacings_Edge_1');
	q{2} = getinfo(reply, 'Spacings_Edge_2');
	q{3} = getinfo(reply, 'Density_Flag');
	psset(self, 'itsDefaultSpacings', q);
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
	update(self, 1)
end

if nargout > 0, theResult = self; end
