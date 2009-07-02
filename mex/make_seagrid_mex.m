% svn $Id$
%
% Create mex files used by seagrid to generate curvilinear grids.  
% 

cd mexinside;

mex mexinside.c;
test_mexinside;
disp ( 'wait for the figure to finish plotting, then hit any key to continue' );
pause
unix('mv mexinside.mex* ../');

cd ../mexrect;
mex -f ./mexopts.sh -v mexrect.F
test_mexrect;
disp ( 'hit any key to continue' );
pause
unix('mv mexrect.mex* ../');

cd ../mexsepeli;
mex -f ./mexopts.sh -v -o mexsepeli mexgateway_sepeli.c sepeli.F
test_mexsepeli;
unix('mv mexsepeli.mex* ../');

cd ../



