function theResult = doevent(self, theCommand, theMessage)

% seagrid/doevent -- Event handler for "seagrid" object.
%  doevent(self, theCommand, theMessage) processes theCommand
%   and theMessage on behalf of self, a "seagrid" object.
 
% svn $Id$
%=======================================================================
% Copyright (C) 1999 Dr. Charles R. Denham, ZYDECO.
%  All Rights Reserved.
%   Disclosure without explicit written consent from the
%    copyright owner does not constitute publication.
%=======================================================================
 
% Version of 07-Apr-1999 20:44:42.
% Updated    23-Apr-2001 16:00:02.

if nargin < 1, help(mfilename), return, end
if nargin < 2, theCommand = ''; end
if nargin < 3, theMessage = ''; end

switch lower(theCommand)
case 'callback'
	if ~isempty(gcbo)
		switch get(gcbo, 'Type')
		case {'uimenu', 'uicontrol'}
			theCommand = get(gcbo, 'Tag');
		end
	end
otherwise
end

theCommand = translate(self, theCommand);

if nargout > 0, theResult = self; end

if verbose(self)
	disp([' ## ' mfilename ' ' theCommand])
end

theFigure = ps(self);

switch lower(theCommand)
case 'buttondownfcn'
	switch get(gcbo, 'Type')
	case 'surface'
		theMaskTool = psget(self, 'itsMaskTool');
		if isequal(get(gcbo, 'Tag'), 'masktool')
			feval('grid2mask', theCommand);
			result = feval('grid2mask', 'mask');
			theMask = isnan(result);
			theMask = ~~theMask;
			theLand = theMask;
			theWater = ~theLand;
			psset(self, 'itsMask', theMask)
			psset(self, 'itsWater', theWater)
			psset(self, 'itsLand', theLand)
			bathy = psget(self, 'itsGriddedBathymetry');
		end
	end
case {'deletedepths', 'setalltozero'}
	theMaskTool = psget(self, 'itsMaskTool');
	if any(theMaskTool)
		domasktool(self)   % Delete it.
	end
	psset(self, 'itsGriddedBathymetry', [])
	psset(self, 'itsBathymetryFlag', 0)
	h = findobj(gcf, 'Type', 'uimenu', 'Tag', 'Grid Depths');
	set(h, 'Checked', 'off')
	doupdate(self, 1)
case 'deletemask'
	theMaskTool = psget(self, 'itsMaskTool');
	if any(theMaskTool)
		domasktool(self)   % Delete it.
	end
	psset(self, 'itsMask', [])
	psset(self, 'itsLand', [])
	psset(self, 'itsWater', [])
	psset(self, 'itsMaskingFlag', 0)
	psset(self, 'itsMaskToolFlag', 0)
	h = findobj(gcf, 'Type', 'uimenu', 'Tag', 'Land Mask');
	set(h, 'Checked', 'off')
	doupdate(self, 1)
case {'depths', 'landmask', 'depthsandlandmask'}
	switch lower(theCommand)
	case 'depths'
		psset(self, 'itsBathymetryFlag', 1)
	case 'landmask'
		psset(self, 'itsMaskingFlag', 1)
	case 'depthsandlandmask'
		psset(self, 'itsBathymetryFlag', 1)
		psset(self, 'itsMaskingFlag', 1)
	end
	psset(self, 'itsMaskToolFlag', 1)
	needsUpdate = 1;
	doupdate(self, needsUpdate)
case 'settoallland'
	theMask = psget(self, 'itsMask');
	if any(any(theMask))
		theMaskTool = psget(self, 'itsMaskTool');
		if any(theMaskTool), domasktool(self), end
		theMask = psget(self, 'itsMask');
		theMask = ~zeros(size(theMask));   % All Land.
		theLand = theMask;
		theWater = ~theLand;
		psset(self, 'itsMask', theMask)
		psset(self, 'itsLand', theLand)
		psset(self, 'itsWater', theWater)
		psset(self, 'itsMaskingFlag', 0)
		h = findobj(theFigure, 'Type', 'uimenu', 'Label', 'Land Mask');
		set(h, 'Checked', 'off')
		doupdate(self, 0)
	end
case 'settoallwater'
	theMask = psget(self, 'itsMask');
	if any(any(theMask))
		theMaskTool = psget(self, 'itsMaskTool');
		if any(theMaskTool), domasktool(self), end
		theMask = psget(self, 'itsMask');
		theMask = ~~zeros(size(theMask));   % All water.
		theLand = theMask;
		theWater = ~theLand;
		psset(self, 'itsMask', theMask)
		psset(self, 'itsLand', theLand)
		psset(self, 'itsWater', theWater)
		psset(self, 'itsMaskingFlag', 0)
		h = findobj(theFigure, 'Type', 'uimenu', 'Label', 'Land Mask');
		set(h, 'Checked', 'off')
		doupdate(self, 0)
	end
case {'corner', 'down', 'motion', 'up'}
	getboundary(self, theCommand, theMessage)
case 'help'
	dohelp(self)
case 'aboutseagrid'
	help_about(self)
case 'bugs'
	help_bugs(self)
case 'depthsdisplay'
	help_depths(self)
case 'converttoecom'
	help_ecom(self)
case 'converttoroms'
	help_roms(self)
%	case 'converttoscrum'
%		help_scrum(self)
case 'futurework'
	help_future(self)
case 'general'
	help_general(self)
case 'gettingstarted'
	help_start(self)
case 'wwwhomepage'
	sgweb
case 'loading'
	help_loading(self)
case 'mathematics'
	help_math(self)
case 'masking'
	help_masking(self)
case 'menus'
	help_menus(self)
case 'orthogonalitydisplay'
	help_orthogonality(self)
case 'points'
	help_points(self)
case 'private'
	help_private(self)
case 'projection'
	help_projection(self)
case 'saving'
	help_saving(self)
case 'setupdialog'
	help_setup_dialog(self)
case 'spacing'
	help_spacing(self)
case 'units'
	help_units(self)
case 'updating'
	help_updating(self)
case 'version'
	version(self)
case 'warranty'
	help_warranty(self)
case 'wwwhomepage'
	web(self)
case 'zooming'
	help_zooming(self)
case 'setup'
	dosetup(self)
	doupdate(self)
case 'new'
	seagrid
case 'coastline'
	h = findobj('Type', 'line', 'Tag', 'coastline');
	showhide(h)
case {'bathymetry', 'bathymetricpoints'}
	h = findobj('Type', 'line', 'Tag', 'bathymetry');
	showhide(h)
case 'depths'
	h = findobj('Type', 'line', 'Tag', 'griddepths');
	showhide(h)
case 'mask'
	h = findobj('Type', 'line', 'Tag', 'landmask');
	showhide(h)
case 'coastlinefile'
	[theFile, thePath] = uigetfile('*.*', 'Select a Coastline File:');
	if ~any(theFile), return, end
	if thePath(end) ~= filesep, thePath = [thePath filesep]; end
	theCoastlineFile = [thePath theFile];
	psset(self, 'itsCoastlineFile', theCoastlineFile)
	getcoastline(self)
	doticks(self)
case 'bathymetryfile'
	[theFile, thePath] = uigetfile('*.*', 'Select a Bathymetry File:');
	if ~any(theFile), return, end
	if thePath(end) ~= filesep, thePath = [thePath filesep]; end
	theBathymetryFile = [thePath theFile];
	psset(self, 'itsBathymetryFile', theBathymetryFile)
	getbathymetry(self)
case 'seagridfile'
	[theFile, thePath] = uigetfile('*.mat', 'Select SeaGrid Input:');
	if ~any(theFile), return, end
	if thePath(end) ~= filesep, thePath = [thePath filesep]; end
	theSeaGridInputFile = [thePath theFile];
	psset(self, 'itsSeaGridInputFile', theSeaGridInputFile)
	doload(self, theSeaGridInputFile)
	thePoints = psget(self, 'itsPoints');
	if length(thePoints) >= 4
		set(gca, 'ButtonDownFcn', '', 'UserData', [])
	end
case 'boundaryfile'
	[theFile, thePath] = uigetfile('*.*', 'Select a Boundary File:');
	if ~any(theFile), return, end
	if thePath(end) ~= filesep, thePath = [thePath filesep]; end
	theBoundaryFile = [thePath theFile];
	psset(self, 'itsBoundaryFile', theBoundaryFile)
	thePoints = load(theBoundaryFile);
% Corrected on 01/14/2000 to project loaded boundary.
% Mike Whitney <whitney@neward.cms.udel.edu>.
	lon = thePoints(:,1); lat = thePoints(:,2);
	theProjection = psget(self, 'itsProjection');
	switch theProjection   % Needs cleaning up.
	case {'none', 'Geographic'}
		theProjection = 'Geographic';
		x = lon; y = lat;
	otherwise
		sg_proj(theProjection)
		[x, y] = sg_ll2xy(lon, lat);
		x = real(x);
		y = real(y);
	end
	thePoints(:,1) = x; thePoints(:,2) = y;
	psset(self, 'itsPoints', thePoints)
	doupdate(self, 1)
case 'boundaryfile_XXX'
	[theFile, thePath] = uigetfile('*.*', 'Select a Boundary File:');
	if ~any(theFile), return, end
	if thePath(end) ~= filesep, thePath = [thePath filesep]; end
	theBoundaryFile = [thePath theFile];
	psset(self, 'itsBoundaryFile', theBoundaryFile)
	thePoints = load(theBoundaryFile);
	psset(self, 'itsPoints', thePoints)
	doupdate(self, 1)
case 'ecomoutput_XXX'
	if (0)
		[theFile, thePath] = uiputfile('ecom.out', 'Save As Ecom File:');
		if ~any(theFile), return, end
		if thePath(end) ~= filesep, thePath = [thePath filesep]; end
		theEcomFile = [thePath theFile];
		psset(self, 'itsEcomFile', theEcomFile)
	end
case 'scrumoutput_XXX'
	if (0)
		[theFile, thePath] = uiputfile('scrum.out', 'Save As Scrum File:');
		if ~any(theFile), return, end
		if thePath(end) ~= filesep, thePath = [thePath filesep]; end
		theScrumFile = [thePath theFile];
		psset(self, 'itsScrumFile', theScrumFile)
	end
case 'save'
	theSeaGridOutputFile = psget(self, 'itsSeaGridOutputFile');
	dosave(self, theSeaGridOutputFile)
case 'saveas'
	dosave(self)
case 'saveasjpeg'
	figure(ps(self))
	fig2jpeg
case 'reverttosaved'
	theSeaGridOutputFile = psget(self, 'itsSeaGridOutputFile');
	doload(self, theSeaGridOutputFile)
	
case 'controlpoints'   % Show control-points.
	h = psget(self, 'itsSpacers');
	set([h{1} h{2}], 'Visible', 'off')
	h = psget(self, 'itsPoints');
	set(h, 'Visible', 'on')
	psset(self, 'itsPointFlag', 1)
	psset(self, 'itsSpacerFlag', 0)
	h = findobj('Type', 'line', 'Tag', 'edge');
%	theButtonDownFcn = 'psevent(ps(theFigure), ''down'')';
	theButtonDownFcn = 'psevent down';
	if any(h), set(h, 'ButtonDownFcn', theButtonDownFcn), end
	
case 'spacers'   % Show spacers.
	h = psget(self, 'itsPoints');
	set(h, 'Visible', 'off')
	h = psget(self, 'itsSpacers');
	h = [h{1}(2:end-1) h{2}(2:end-1)];   % Ends not visible.
	set(h, 'Visible', 'on')
	psset(self, 'itsPointFlag', 0)
	psset(self, 'itsSpacerFlag', 1)
	h = findobj('Type', 'line', 'Tag', 'edge');
	theButtonDownFcn = '';
	if any(h), set(h, 'ButtonDownFcn', theButtonDownFcn), end
case 'spacercount'   % Enter spacer counts.
	getspacers(self)
case 'spacingsetup'   % Enter spacing functions.
	getspacings(self)
	
case 'orthogonality'
	theOrthogonalityFlag = psget(self, 'itsOrthogonalityFlag');
	if isempty(theOrthogonalityFlag), theOrthogonalityFlag = 0; end
	switch theOrthogonalityFlag
	case 0
		theOrthogonalityFlag = 1;
		set(gcbo, 'Checked', 'on')
	otherwise
		theOrthogonalityFlag = 0;
		set(gcbo, 'Checked', 'off')
	end
	psset(self, 'itsOrthogonalityFlag', theOrthogonalityFlag)
	doupdate(self)
case 'coriolis-XXX'
	disp([' ## ' theCommand ' ' theMessage])
case 'roll'
	doroll(self)
case 'axesequal'
	axis equal
case 'margin'
	zoomsafe(0.9)
case 'graticule'
	switch psget(self, 'itsGraticule')
	case 'on'
		psset(self, 'itsGraticule', 'off')
		set(gcbo, 'Checked', 'off')
	case 'off'
		psset(self, 'itsGraticule', 'on')
		set(gcbo, 'Checked', 'on')
	otherwise
	end
	dograticule(self)
case 'menubar'
	switch get(theFigure, 'MenuBar')
	case 'figure'
		set(theFigure, 'MenuBar', 'none')
	otherwise
		set(theFigure, 'MenuBar', 'figure')
	end
case 'zoomin'
	dograticule(self, 'off')
	zoomsafe(2)
	doticks(self)
case 'zoomout'
	dograticule(self, 'off')
	zoomsafe(0.5)
	doticks(self)
case 'nozoom'
	dograticule(self, 'off')
	zoomsafe out
	doticks(self)
case 'flat'
	shading flat
	refresh
case 'interp'
	shading interp
case 'print'
	printsafe(gcf, '-v')
case 'degrees'
	psset(self, 'itsMapUnits', 'degrees')
	doticks(self)
case 'kilometers'
	psset(self, 'itsMapUnits', 'kilometers')
	doticks(self)
case 'projected'
	psset(self, 'itsMapUnits', 'projected')
	doticks(self)
case 'showdetails'
	assignin('base', 'ans', self)
	evalin('base', 'ans')
case 'verbose'
	switch get(gcbo, 'Checked')
	case 'on'
		verbose(self, ~~0)
		set(gcbo, 'Checked', 'off')
	otherwise
		verbose(self, ~~1)
		set(gcbo, 'Checked', 'on')
	end
case 'refreshscreen'
	refresh(ps(self))
case 'update'
	doupdate(self, 1)
case 'resizefcn'
	doticks(self)
case {'closerequestfcn', 'quit'}
	doquit(self)
otherwise
	doevent(super(self), theCommand, theMessage)   % Very important.
end
