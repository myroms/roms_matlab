function help_menus(self)

% SeaGrid Menu Layout
%    (Accelerator keys shown in parentheses)
%
% <SeaGrid>
%    About SeaGrid...        	Purpose and authorship.
%    New...                  	New SeaGrid window.
%    Load...                 	Load data files.
%       Coastline File...       Load (lon, lat) coastline data.
%       Bathymetry File...      Load (lon, lat, depth) data.
%       Boundary File...        Load (lon, lat, isCorner) data.
%       SeaGrid File...         Load grid from "SeaGrid" file.
%    Save                    	Save grid to current SeaGrid file.
%    Save As...              	Save grid to new SeaGrid file.
%    Revert To Saved         	Revert to last saved SeaGrid file.
%    Print...                	Print SeaGrid figure.
%    Quit                    	Quit SeaGrid window.
%
% <View>
%    Setup... (')               Dialog for SeaGrid parameters.
%    Control Points          	Show draggable boundary points.
%    Roll                    	Roll principle corner CCW.
%    Spacers                 	Show draggers for grid spacing.
%    Spacer Count...         	Adjust number of spacers.
%    Spacing Setup...        	Set default spacing algorithm.
%    Axes Equal (=)             Aspect ratio = 1.
%    Margin (-)                 Make small margin.
%    Zoom In ([)                Zoom map in x2.
%    Zoom Out (])               Zoom map out x2.
%    No Zoom (\)                Size map to fill figure.
%    Map Units               	Set map label units.
%       Degrees                 Degrees.
%       Kilometers              Kilometers.
%       Projected               Raw projected units.
%    Show/Hide               	Toggle visibility.
%       Coastline               Toggle coastline.
%       Bathymetry              Toggle bathymetry points.
%       Graticule (G)           Toggle map graticule.
%       MenuBar (M)             Toggle Matlab menubar.
%    Details                 	Display SeaGrid properties.
%    Refresh                    Redraw without recomputation.
%    Update (U)                 Full computational update.
%
% <Compute>
%    Depths And Auto-Mask		Depth and automatic land-mask.
%    Depths                     Depths only.
%    Auto-Mask                  Automatic land-mask only.
%    Erase Mask              	Clear land mask.
%    Orthogonality           	Show orthogonality error.
%
% <Help>
%    Help...                    Dialog for scrollable help.
%    Bugs                       Known bugs and features.
%    Convert To Ecom...         SeaGrid format --> ECOM format.
%    Convert To Scrum...        SeaGrid format --> SCRUM format.
%    Depths Display...          Bathymetry gridding and masking.
%    Future Work...             To do.
%    General...                 General remarks.
%    Getting Started...         Easy startup.
%    Loading...                 Load data and grids.
%    Mathematics...             Algorithm.
%    Masking...                 Land masking.
%    Menus...                   Menu layout.
%    Orthogonality Display...   Angular errors in grid.
%    Points...                  Manipulate control points.
%    Private...                 Private functions.
%    Projection...              Select map projection.
%    Save...                    Save grid to file.
%    Setup Dialog...            Dialog for major parameters.
%    Spacing...                 Manipulate grid spacing.
%    Units...                   Map units.
%    Updating...                Fully update the grid.
%    Version...                 Version of this program.
%    Zooming...                 Zoom via mouse clicks.
%    Warranty...                Caveat emptor.
%    WWW Home Page...           Go to SeaGrid WWW Hole Page.

% svn $Id$
%=======================================================================
% Copyright (C) 2000 Dr. Charles R. Denham, ZYDECO.
%  All Rights Reserved.
%   Disclosure without explicit written consent from the
%    copyright owner does not constitute publication.
%=======================================================================
 
% Version of 04-Jan-2000 09:26:06.
% Updated    17-Aug-2000 10:53:59.

hint(help(mfilename), 'SeaGrid Menus')
