% Simple test of mexsepelli
%

% svn $Id$
%=======================================================================
% Copyright (C) 1997 Dr. Charles R. Denham, ZYDECO.
%  All Rights Reserved.
%   Disclosure without explicit written consent from the
%    copyright owner does not constitute publication.
%=======================================================================

load mexsepeli_testdata;

[xgrid,ygrid] = mexsepeli ( x, y, l2, m2, seta, sxi );


xdiff = xgrid(:) - goodx(:);
ydiff = ygrid(:) - goody(:);

figure;
subplot(2,1,1);
plot(xdiff); 
ylabel ( 'x grid error' );
title('Test of Mexsepeli:  difference between computed grid and solution grid, should be all zero' );
subplot(2,1,2);
plot ( ydiff );
ylabel ( 'y grid error' );
xlabel ( 'grid index' );
