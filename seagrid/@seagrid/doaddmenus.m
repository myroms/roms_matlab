function theResult = doaddmenus(self)

% seagrid/doaddmenus -- Add menus to "seagrid" window.
%  doaddmenus(self) adds menus to the window associated
%   with self, a "seagrid" object.
%
% (Now Presto Compliant)
 
% svn $Id$
%=======================================================================
% Copyright (C) 1997 Dr. Charles R. Denham, ZYDECO.
%  All Rights Reserved.
%   Disclosure without explicit written consent from the
%    copyright owner does not constitute publication.
%=======================================================================
 
% Version of 21-Jul-1997 16:30:13.
% Updated    23-Apr-2001 12:03:44.

if nargin < 1, help(mfilename), return, end

theFigure = ps(self);

theMenuBar = get(theFigure, 'MenuBar');

h = findobj(gcf, 'Type', 'uimenu');
if any(h), delete(h), end

% New menus.

theLabels = ...
{ ...
'<SeaGrid>'
'>About SeaGrid...'
'-New'
'>Load...'
'>>Coastline File...'
'>>Bathymetry File...'
'->Boundary File...'
'->SeaGrid File...'
'-Save'
'>Save As...'
'>Revert To Saved'
'-Save As JPEG...'
'-Print...'
'-Quit'
...
'<View>'
{'>Setup...', ''''}
'-Control Points'
'>Roll'
'-Spacers'
'>Spacer Count...'
'>Spacing Setup...'
{'-Axes Equal', '='}
{'>Margin', '-'}
{'-Zoom In', '['}
{'>Zoom Out', ']'}
{'>No Zoom', '\'}
'-Map Units'
'>>Degrees'
'>>Kilometers'
'->Projected'
...
'<Compute>'
'>Depths And Land Mask'
'-Depths'
% '>Grid Depths Once'
% '>Grid Depths Always'
'>Delete Depths'
'-Land Mask'
% '>Land Mask Once'
% '>Land Mask Always'
'>Set To All Land'
'>Set To All Water'
'>Delete Mask'
%	'>Erase Mask'
'-Orthogonality'
{'-Update', 'U'}
...
'<Toggle>'
'>Coastline'
'>Bathymetric Points'
%	'-Depths'
%	'>Mask'
{'-Graticule', 'G'}
{'-MenuBar', 'M'}
'-Verbose'
'-Show Details'
'-Refresh Screen'
...
'<Help>'
{'>Help...', 'H'}
'-Bugs...'
'>Depths Display...'
'>Convert To Ecom...'
'>Convert To Roms...'
'>Future Work...'
'>General...'
'>Getting Started...'
'>Loading...'
'>Masking...'
'>Mathematics...'
'>Menus...'
'>Orthogonality Display...'
'>Points...'
'>Private...'
'>Projection...'
'>Saving...'
'>Setup Dialog...'
'>Spacing...'
'>Units...'
'>Updating...'
'>Version...'
'>Zooming...'
'-Warranty...'
'>WWW Home Page...'
};

theMenus = menu(self, theLabels);

theEventHandlers = { ...
	'aboutseagrid', 'help_about', ...
	'new', 'donew', ...
	'coastlinefile', 'getcoastline', ...
	'bathymetryfile', 'getcoastline', ...
	'boundaryfile', 'getboundary', ...
	'seagridfile', 'getseagrid', ...
	'save', 'dosave', ...
	'saveas', 'dosaveas', ...
	'reverttosaved', 'dorevert', ...
	'print', 'doprint', ...
	'quit', 'doquit', ...
	
};

self = handler(self, theEventHandlers{:});

enable(self)

set(theFigure, 'MenuBar', 'none')

if nargout > 0, theResult = theMenus; end
