function theResult = getspacers(self)

% seagrid/getspacers -- Get spacer counts.
%  getspacers(self) invokes a dialog to get the
%   number of spacers along edges #1 and #2.
 
% svn $Id$
%=======================================================================
% Copyright (C) 1999 Dr. Charles R. Denham, ZYDECO.
%  All Rights Reserved.
%   Disclosure without explicit written consent from the
%    copyright owner does not constitute publication.
%=======================================================================
 
% Version of 03-Aug-1999 10:37:40.
% Updated    10-Jan-2000 16:33:03.

if nargin < 1, help(mfilename), return, end

theGridSize = psget(self, 'itsGridSize');

theSpacerCounts = psget(self, 'itsSpacerCounts');

if isempty(theSpacerCounts)
	theSpacerCounts = [5 5];
end

d.Spacers_Edge_1 = theSpacerCounts(1);   % uigetinfo dialog structure.
d.Spacers_Edge_2 = theSpacerCounts(2);
SeaGrid_Spacer_Setup = d;

if (0)
%	reply = uigetinfo(SeaGrid_Spacer_Setup);
else
	reply = guido(SeaGrid_Spacer_Setup);
end

if ~isempty(reply)
	theSpacerCounts(1) = getinfo(reply, 'Spacers_Edge_1');
	theSpacerCounts(2) = getinfo(reply, 'Spacers_Edge_2');
	theSpacerCounts = max(theSpacerCounts, 1);
	psset(self, 'itsSpacerCounts', theSpacerCounts)
	doupdate(self, 1)
end

if nargout > 0, theResult = self; end
