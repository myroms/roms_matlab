function theResult = doreproject(self)

% seagrid/doreproject -- Re-project the SeaGrid map.
%  doreproject(self) re-projects the SeaGrid map from
%   the current projection to the new projection of
%   record.
 
% svn $Id$
%=======================================================================
% Copyright (C) 1999 Dr. Charles R. Denham, ZYDECO.
%  All Rights Reserved.
%   Disclosure without explicit written consent from the
%    copyright owner does not constitute publication.
%=======================================================================
 
% Version of 28-Dec-1999 17:03:30.
% Updated    28-Dec-1999 17:03:30.

if nargout > 0, theResult = self; end

theProjection = psget(self, 'itsProjection');
theProjectionCenter = psget(self, 'itsProjectionCenter');

theNewProjection = psget(self, 'itsNewProjection');
theNewProjectionCenter = psget(self, 'itsNewProjectionCenter');

if isequal(theProjection, theNewProjection) & ...
	isequal(theProjectionCenter, theNewProjectionCenter)
	return
end

theFigure = ps(self);
figure(theFigure)
h = get(gca, 'Children');

for i = 1:length(h)
	switch get(h(i), 'Type')
	case {'line', 'patch', 'surface'}
		x = get(h(i), 'XData');
		y = get(h(i), 'YData');
		sg_proj(theProjection, theProjectionCenter)
		[lon, lat] = sg_xy2ll(x, y);
		sg_proj(theNewProjection, theNewProjectionCenter)
		[x, y] = sg_ll2xy(lon, lat);
		set(h(i), 'XData', x, 'YData', y)
	otherwise
	end
end

psset(self, 'itsProjection', theNewProjection)
psset(self, 'itsProjectionCenter', theNewProjectionCenter)

psset(self, 'itsOldProjection', [])
psset(self, 'itsOldProjectionCenter', [])

if nargout > 0, theResult = self; end
