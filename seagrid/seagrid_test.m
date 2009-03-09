function seagrid_demo

% seagrid_demo -- Demo of "seagrid" facility.
 
% svn $Id$
%=======================================================================
% Copyright (C) 1999 Dr. Charles R. Denham, ZYDECO.
%  All Rights Reserved.
%   Disclosure without explicit written consent from the
%    copyright owner does not constitute publication.
%=======================================================================
 
% Version of 28-Apr-1999 13:10:18.
% Updated    24-Jul-2000 13:26:43.

clear all

sg

cd test_data

theCoastline = which('amazon_coast.mat')
theBathymetry = which('amazon_bathy.mat')

cd ..

self = seagrid(theCoastline, theBathymetry)

assignin('base', 'self', self)
