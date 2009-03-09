function seagrid_bundle

% seagrid_bundle -- Bundler for "SeaGrid" toolbox.
%  seagrid_bundle (no argument) bundles the files
%   that comprise the "SeaGrid" toolbox.  The
%   result is "seagrid_install.p".
 
% svn $Id$
%=======================================================================
% Copyright (C) 1999 Dr. Charles R. Denham, ZYDECO.
%  All Rights Reserved.
%   Disclosure without explicit written consent from the
%    copyright owner does not constitute publication.
%=======================================================================
 
% Version of 15-Apr-1999 16:24:49.
% Updated    25-Oct-2000 15:33:34.

if (0)
	ELAPSED = 'disp(['' ## Elapsed time: '' int2str(toc) '' s''])   % Timing.';
	tic
end

theDirs = {
	'seagrid'
	'seagrid:@seagrid'
	'seagrid:test_data'
};

theMFiles = {

	mfilename
	
	'README'
	
	'seagrid_test'
	'grid_test'
	'sg'
	'hello'
	'setdef'
	
	'sg_mercator'
	'sg_stereographic'
	'sg_lambert_equal_area'
	
	'sg_proj'
	'sg_ll2xy'
	'sg_xy2ll'

	'insidesafe'
	'labelsafe'
	'printsafe'
	'setsafe'
	'splinesafe'
	'zoomsafe'

	'fieldrename'
	'pdf2tick'

	'linticks'

	'showhide'

	'busy'
	'idle'

	'ismac'
	'ispc'
	'vi'

	'rescale'
	'respace'

	'earthdist'
	'splinep'
	'splineq'

	'fpt'
	'fps'
	'rect'
	'rect2grid'
	'mexrect2grid'

	'grid2mask'
	
	'griddata1'
	'tgriddata'

	'monocline'

	'seagrid2ecom'
	'seagrid2roms'
	'seagrid2scrum'

	'seagrid_helpdlg'

	'fig2jpeg'
	'guido'
	'hint'
	'getinfo'
	'setinfo'
};

theMFiles = sort(theMFiles);

okay = 1;
for i = 1:length(theMFiles)
	if isempty(which(theMFiles{i}))
		disp([' ## Not found: "' theMFiles{i} '"'])
		okay = 0;
	end
end

if ~okay
	disp(' ## Unable to continue.  Some files are missing.')
	return
end

theClasses = {
	'seagrid'
};

theMatFiles = {
	'amazon_coast.mat'
	'amazon_bathy.mat'
};

theMessages = {
	' '
	' ##  Requires the "presto" toolbox (Denham)'
	' '
	' ##  Unix and PC require these Mex-files:'
	' ##   from "gridgen" (Evans):'
	' ##       "mexinside"'
	' ##       "mexrect"'
	' ##       "mexsepeli"'

	' '
	' ## To get started, put the "seagrid" and "presto" folders'
	' ## (and other relevant files) in your Matlab path.'
	' ## Restart Matlab, then execute "seagrid" at the Matlab prompt.'
};

if (0)

fclose('all');

setdef(mfilename)

newversion seagrid

dst = 'seagrid_install';

self = bundle(dst, 'new');

disp(' ')
disp(self)

self = add_checkdir(self, theDirs);

self = add_setdir(self, 'seagrid');

self = add_mfile(self, theMFiles);
self = add_class(self, theClasses);

% Add test_data.

self = add_setdir(self, 'test_data');
cd test_data
self = add_binary(self, theMatFiles);
self = add_setdir(self, '..');
cd ..

self = add_setdir(self, '..');

self = add_message(self, theMessages);

self = make_pcode(self);

eval(ELAPSED)

end

% New bundling scheme, using BUND.

at(mfilename)

newversion seagrid

bund new seagrid
bund setdir seagrid
bund('mfile', theMFiles)
bund('class', theClasses)
cd test_data
bund setdir test_data
bund('binary', theMatFiles)
cd ..
bund cd ..
bund cd ..
bund('message', theMessages)
bund close
