function theResult = seagrid(theCoastlineFile, theBathymetryFile, theProjection)

% seagrid/seagrid -- Constructor for "seagrid" class.
%  seagrid('demo') shows itself using the Amazon "test_data".
%  seagrid('theCoastlineFile', 'theBathymetryFile', 'theProjection')
%   constructs a "seagrid" interface, using the given coastline file,
%   bathymetry file, and projection (default = 'Mercator').
%  seagrid('version') reports the current SeaGrid version number.
 
% svn $Id$
%=======================================================================
% Copyright (C) 1999 Dr. Charles R. Denham, ZYDECO.
%  All Rights Reserved.
%   Disclosure without explicit written consent from the
%    copyright owner does not constitute publication.
%=======================================================================
 
% Version of 07-Apr-1999 15:44:33.
% Updated    13-Feb-2001 10:57:17.

if nargin == 1 & isequal(theCoastlineFile, 'demo')
	oldPWD = pwd;
	setdef(mfilename)
	cd ..
	cd test_data
	try
		result = seagrid('amazon_coast.mat', 'amazon_bathy.mat');
	catch
		result = [];
	end
	cd(oldPWD)
	if nargout > 0
		theResult = result;
	else
		assignin('base', 'ans', result)
	end
	return
end

if nargout > 0, theResult = []; end
if nargin < 1, theCoastlineFile = ''; end
if nargin < 2, theBathymetryFile = ''; end
if nargin < 3, theProjection = 'Mercator'; end

% Locate mex-files for Unix and PCWIN.

c = lower(computer);
isMac = isequal(c(1:3), 'mac');
if ~isMac
	hasMex = (exist('mexrect', 'file') == 3) & ...
				(exist('mexsepeli', 'file') == 3) & ...
				(exist('mexinside', 'file') == 3);
	if ~hasMex
		disp(' ')
		disp(' ## Unable to locate Mex-files.  Please install')
		disp(' ## "mexrect", "mexsepeli", and "mexinside",')
		disp(' ## then adjust your Matlab path accordingly.')
		disp(' ')
	end
end

% Locate libraries.

okay = 1;
if ~any(which('ps/ps'))
	disp(' ## SeaGrid requires the "Presto" toolbox.')
	okay = 0;
end

if 0 & ~any(which('m_proj'))   % No longer used.
	disp(' ## SeaGrid requires the "M_Map" toolbox.')
	okay = 0;
end

if ~okay
	disp(' ## Please install and/or include in Matlab path.')
	return
end

% Version.

if isequal(lower(theCoastlineFile), 'version')
	help seagrid/version
	return
end

% Is PWD in a class-directory.

%if (isvms & any(pwd == '#')) | (~isvms & any(pwd == '@'))
%	disp([' ## "' pwd '" is a class-directory.'])
%	disp([' ## Please start SeaGrid from outside any class-directory.'])
%	return
%end

if (any(pwd == '@'))
	disp([' ## "' pwd '" is a class-directory.'])
	disp([' ## Please start SeaGrid from outside any class-directory.'])
	return
end

busy

theFigure = figure('Name', 'SeaGrid');

theStruct.ignore = [];
self = class(theStruct, 'seagrid', ps(theFigure));
psbind(self)

psset(self, 'itsCoastlineFile', theCoastlineFile)
psset(self, 'itsBathymetryFile', theBathymetryFile)

theProjectionCenter = [0 0 0];
theLongitudeBounds = [];
theLatitudeBounds = [];
theMapUnits = 'degrees';

theFigure = ps(self);

theMenus = doaddmenus(self);

theInterpFcn = 'splinesafe';
theInterpMethod = 1;
theGriddingMethod = 'linear';   % For bathymetry.

self = doinitialize(self, ...
					'itNeedsUpdate', ~~0, ...
					'itIsVerbose', ~~0, ...
					'itsCoastlineFile', theCoastlineFile, ...
					'itsCoastlineColor', [4 2 1]/4, ...
					'itsBathymetryFile', theBathymetryFile, ...
					'itsBathymetryColor', [4 4 5]/5, ...
					'itsProjection', theProjection, ...
					'itsProjectionCenter', theProjectionCenter, ...
					'itsLongitudeBounds', theLongitudeBounds, ...
					'itsLatitudeBounds', theLatitudeBounds, ...
					'itsMapUnits', theMapUnits, ...
					'itsGraticule', 'off', ...
					'itsGridSize', [10 10], ...
					'itsEraseMode', 'xor', ...
					'itsCornerColor', 'r', ...
					'itsCornerMarker', 'o', ...
					'itsCornerLineStyle', 'none', ...
					'itsCornerTag', 'corner-point', ...
					'itsEdgeColor', 'b', ...
					'itsEdgeMarker', 'o', ...
					'itsEdgeLineStyle', '-', ...
					'itsEdgePointTag', 'edge-point', ...
					'itsEdgeTag', 'edge', ...
					'itsGriddedBathymetry', [], ...
					'itsMask', [], ...
					'itsLand', [], ...
					'itsWater', [], ...
					'itsOrthogonalityTag', 'orthogonality', ...
					'itsGridLineColor', 'g', ...
					'itsEndSlopeFlag', 1, ...
					'itsClippingDepths', [0 Inf], ...
					'itsPointFlag', 1, ...
					'itsSpacerFlag', 0, ...
					'itsSpacerMarker', 's', ...
					'itsSpacerCounts', [5 5], ...
					'itsSpacedEdges', [1 2], ...
					'itsDefaultSpacings', {'s', 's'}, ...
					'itsInterpFcn', theInterpFcn, ...
					'itsInterpMethod', theInterpMethod, ...
					'itsGriddingMethod', theGriddingMethod, ...
					'itUsesMexFile', ~~1);

if nargin > 0
	load(theCoastlineFile)
	if exist('s', 'var') == 1 & isstruct(s)
		theSeagridFile = theCoastlineFile;
		doload(self, theSeagridFile)
		return
	else
		getcoastline(self)
	end
end

if nargin > 1, getbathymetry(self), end

enable(self)
set(theFigure, 'WindowButtonDownFcn', '')

f = findobj(gcf, 'Type', 'uimenu', 'Parent', gcf);
if any(f), set(f, 'Callback', ''), end

doticks(self)

getboundary(self)

if nargout > 0, theResult = self; end
