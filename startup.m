function startup

% startup -- User script configuration for Matlab.  It can set default
%            paths, define Handle Graphics defaults, or predefine
%            variables in your workspace.

% svn $Id$
%===========================================================================%
%  Copyright (c) 2002-2012 The ROMS/TOMS Group                              %
%    Licensed under a MIT/X style license                                   %
%    See License_ROMS.txt                           Hernan G. Arango        %
%===========================================================================%

global IPRINT
IPRINT=0;

format long g

% Change "my_root" to the appropriate path were these matlab scripts are
% installed in your computer.

my_root='~/ocean/repository';

path(path, fullfile(my_root, 'matlab', '4dvar', ''))
path(path, fullfile(my_root, 'matlab', 'bathymetry', ''))
path(path, fullfile(my_root, 'matlab', 'boundary', ''))
path(path, fullfile(my_root, 'matlab', 'coastlines', ''))
path(path, fullfile(my_root, 'matlab', 'forcing', ''))
path(path, fullfile(my_root, 'matlab', 'grid', ''))
path(path, fullfile(my_root, 'matlab', 'initial', ''))
path(path, fullfile(my_root, 'matlab', 'landmask', ''))
path(path, fullfile(my_root, 'matlab', 'mex', ''))
path(path, fullfile(my_root, 'matlab', 'netcdf', ''))
path(path, fullfile(my_root, 'matlab', 'seagrid', ''))
path(path, fullfile(my_root, 'matlab', 'seagrid', 'presto', ''))
path(path, fullfile(my_root, 'matlab', 'seawater', ''))
path(path, fullfile(my_root, 'matlab', 't_tide', ''))
path(path, fullfile(my_root, 'matlab', 'tidal_ellipse', ''))
path(path, fullfile(my_root, 'matlab', 'utility', ''))
