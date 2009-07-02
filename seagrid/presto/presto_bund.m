function presto_bund

% presto_bund -- Bundle "presto" software.
%  presto_bund (no argument) bundles the
%   "presto" software into "presto_install.p".
 
% svn $Id$
%=======================================================================
% Copyright (C) 1999 Dr. Charles R. Denham, ZYDECO.
%  All Rights Reserved.
%   Disclosure without explicit written consent from the
%    copyright owner does not constitute publication.
%=======================================================================
 
% Version of 03-Nov-1999 08:11:48.
% Updated    19-Jul-2001 16:46:18.

setdef(mfilename)

theMFiles = {
	mfilename
	'README'
	'setdef'
	'inherit'
	'super'
	'isps'
	'psbind'
	'psevent'
	'ps_test'
};

theMFiles = sort(theMFiles);

theClasses = {
	'presto'
	'ps'
	'pst'
};

theMessage = {
	' '
	' ## To get started, put the "presto" folder in your Matlab'
	' ##  path, then execute "presto_test" at the Matlab prompt.'
};

for i = 1:length(theClasses)
	newversion(theClasses{i})
end

tic

bund new presto
bund setdir presto
bund('mfile', theMFiles)
bund('class', theClasses)
bund cd ..
bund('message', theMessage)
bund close

elapsed
