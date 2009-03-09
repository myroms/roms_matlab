function seagrid2ecom(theSeagridFile, theEcomFile, theSigmaLevels)

% seagrid2ecom -- Convert "seagrid" output to "ecom" input".
%  seagrid2ecom('theSeagridFile', 'theEcomFile') converts
%   the given SeaGrid file to an "ecom" grid file.  Where
%   filenames are not provided, the "uigetfile/uiputfile"
%   dialogs are invoked.  The data represent a grid with
%   twice as many grid lines as needed for ECOM, so that
%   we can work with cell-centers.
%  seagrid(..., ..., theSigmaLevels) uses the given number
%   of levels or vector of levels.
%   The default is (0:-0.1:-1).
 
% svn $Id$
%=======================================================================
% Copyright (C) 1999 Dr. Charles R. Denham, ZYDECO.
%  All Rights Reserved.
%   Disclosure without explicit written consent from the
%    copyright owner does not constitute publication.
%=======================================================================

% Version of 07-Jun-1999 15:43:58.
% Updated    13-Jul-2000 15:43:09.

% Constants.

RCF = 180 / pi;   % Radian conversion factor.
LAND_CODE = -99999;

if nargin < 1, theSeagridFile = '*.mat'; end
if nargin < 2, theEcomFile = 'ecom_model_grid.txt'; end
if nargin < 3, theSigmaLevels = 11; end

if length(theSigmaLevels) == 1 & theSigmaLevels > 1
	theSigmaLevels = linspace(0, 1, theSigmaLevels);
end

% Get filenames.

if any(theSeagridFile == '*')
	[f, p] = uigetfile(theSeagridFile, 'Select SeaGrid File:');
	if ~any(f), return, end
	if p(end) ~= filesep, p(end+1) = filesep; end
	theSeagridFile = [p f]
end

if nargin < 2
	[f, p] = uiputfile(theEcomFile, 'Save to Ecom File:');
	if ~any(f), return, end
	if p(end) ~= filesep, p(end+1) = filesep; end
	theEcomFile = [p f]
end

disp([' ## SeaGrid Source File  : ' theSeagridFile])
disp([' ## Ecom Destination File: ' theEcomFile])

% Load input file.

try
	theSeagridData = load(theSeagridFile, 's');
catch
	disp([' ## Unable to load: "' theSeagridFile '"'])
	return
end

% Open output file.

fecom = fopen(theEcomFile, 'w');

if fecom < 0
	disp([' ## Unable to open: "' theEcomFile '"'])
	disp(' ## If already open elsewhere, please close, then try again.')
	return
end

% How rigid is the following format?
% For large grid cells, the H1 and H2 fields collide.

theFormat = ...
	'%4.0f %4.0f %10.2f %10.2f %10.2f %8.1f %8.1f %8.1f %15.6f %15.6f\n';

s = theSeagridData.s;

lon = s.geographic_grids{1};
lat = s.geographic_grids{2};

% Double the grid-size by simle interpolation.

lon = interp2(lon, 1);
lat = interp2(lat, 1);

% Cell sizes.

if (1)
	dx = earthdist(lon(2:2:end, 1:2:end-1), lat(2:2:end, 1:2:end-1), lon(2:2:end, 3:2:end), lat(2:2:end, 3:2:end));
	dy = earthdist(lon(1:2:end-1, 2:2:end), lat(1:2:end-1, 2:2:end), lon(3:2:end, 2:2:end), lat(3:2:end, 2:2:end));
else
	dx = earthdist(lon(:, 1:end-1), lat(:, 1:end-1), lon(:, 2:end), lat(:, 2:end));
	dy = earthdist(lon(1:end-1, :), lat(1:end-1, :), lon(2:end, :), lat(2:end, :));
end

h1 = dx;   % ECOM nomenclature.
h2 = dy;

% Cell orientations: need to check indices and directions here.

dlon = diff(lon.').';
dlat = diff(lat.').';
clat = cos(lat / RCF);
clat(:, end) = [];
ang = atan2(dlat, dlon .* clat) * RCF;
ang = 0.5 * (ang(2:2:end, 1:2:end) + ang(2:2:end, 2:2:end));

% Geographic cell-centers.

lon = lon(2:2:end-1, 2:2:end-1);
lat = lat(2:2:end-1, 2:2:end-1);

% Bathymetry and mask.

bathymetry = s.gridded_bathymetry;
if isempty(bathymetry)
	bathymetry = zeros(size(lon));
end

mask = s.mask;
if isempty(mask)
	mask = zeros(size(bathymetry));   % All water.
end

% Masks.

mask = ~~mask;   % 1 = land; 0 = water.
land = mask;
water = ~land;

% Impose LAND_CODE around the edge.

land([1 end], :) = LAND_CODE;
land(:, [1 end]) = LAND_CODE;

% Clip and mask the depths.

min_depth = s.clipping_depths(1);
max_depth = s.clipping_depths(2);
bathymetry(bathymetry < min_depth) = min_depth;
bathymetry(bathymetry > max_depth) = max_depth;
bathymetry(land) = LAND_CODE;

% Adjust sizes.

% Which points actually have depths -- centers?
%  Why not keep them all?

if all(size(bathymetry) > size(lon)) & 0   % No longer.
	bathymetry = bathymetry(2:2:end-1, 2:2:end-1);
end

% Which points are actually masked -- centers?
%  Why not keep the entire mask?

if all(size(mask) > size(lon)) & 0   % No longer.
	mask = mask(2:2:end-1, 2:2:end-1);
end

% Coriolis parameter (use latitude).

coriolis = lat;

% Write the number of grid cells.

[rows, columns] = size(lon);

try
	fprintf(fecom, 'Created by SeaGrid on %s\n', datestr(now));
	if (1)
		fprintf(fecom, 'Vertical Segmentation - Sigma Levels (KB)\n');
		sigma = theSigmaLevels(:);   % Column vector;
		sigma = sigma(sigma <= 0 & sigma >= -1);
		if ~any(sigma == 0), sigma = [0; sigma]; end
		if ~any(sigma == -1), sigma = [sigma; -1]; end
		sigma = -sort(-sigma);
		fprintf(fecom, '%s\n', int2str(length(sigma)));
		for i = 1:length(sigma)
			fprintf(fecom, '%4.2f\n', sigma(i));
		end
		fprintf(fecom, ...
			'Horizontal Segmentation (IM & JM) -- Grid Information\n');
	end
	fprintf(fecom, '%5i%5i\n', columns, rows);   % Is the order correct?
catch
	disp([' ## ' lasterr])
	fclose(fecom);
	return
end

% Add a dummy border of masked cells.

if (0)   % OBSOLETE.
	rows = rows + 2;
	columns = columns + 2;
	
	z = zeros(size(lon)+2);
	
	temp = z; temp(2:end-1, 2:end-1) = h1; h1 = temp;
	temp = z; temp(2:end-1, 2:end-1) = h2; h2 = temp;
	temp = z + LAND_CODE; temp(2:end-1, 2:end-1) = bathymetry; bathymetry = temp;
	temp = z; temp(2:end-1, 2:end-1) = ang; ang = temp;
	temp = z; temp(2:end-1, 2:end-1) = coriolis; coriolis = temp;
	temp = z; temp(2:end-1, 2:end-1) = lon; lon = temp;
	temp = z; temp(2:end-1, 2:end-1) = lat; lat = temp;
end

% Output the grid.

for column = 1:columns   % Start at the left.
	remaining = column;
	if remaining == columns | (rem(remaining, 100) == 0)
		disp([' ## Columns remaining: ' int2str(remaining)])
	end
	for row = rows:-1:1   % Start at the bottom.
		try
			i = column;   % ECOM notation: i -> x; j -> y.
			j = rows-row+1;
			fprintf(fecom, theFormat, ...
	          	i, j, ...
			 	h1(row, column), h2(row, column), ...
			  	bathymetry(row, column), ...
			  	ang(row, column), ...
			  	coriolis(row, column), ...
			  	0.0, ...
			  	lon(row, column), lat(row, column));
		catch
			disp([' ## ' lasterr])
			fclose(fecom);
			return
	  	end
	end
end

fclose(fecom);

disp(' ## Done.')
