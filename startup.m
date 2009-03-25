function startup

% startup -- User script configuration for Matlab.  It can set default
%            paths, define Handle Graphics defaults, or predefine
%            variables in your workspace.

% svn $Id$
%===========================================================================%
%  Copyright (c) 2002-2009 The ROMS/TOMS Group                              %
%    Licensed under a MIT/X style license                                   %
%    See License_ROMS.txt                           Hernan G. Arango        %
%===========================================================================%

global IPRINT
IPRINT=0;

format long g

my_root='/home/arango/ocean/toms/repository';

path(path, fullfile(my_root, 'matlab', 'bathymetry', ''))
path(path, fullfile(my_root, 'matlab', 'coastlines', ''))
path(path, fullfile(my_root, 'matlab', 'initial', ''))
path(path, fullfile(my_root, 'matlab', 'landmask', ''))
path(path, fullfile(my_root, 'matlab', 'seagrid', ''))
path(path, fullfile(my_root, 'matlab', 'seawater', ''))
path(path, fullfile(my_root, 'matlab', 'utility', ''))
