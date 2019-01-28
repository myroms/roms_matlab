function startup

% startup -- User script configuration for Matlab.  It can set default
%            paths, define Handle Graphics defaults, or predefine
%            variables in your workspace.

% svn $Id$
%===========================================================================%
%  Copyright (c) 2002-2019 The ROMS/TOMS Group                              %
%    Licensed under a MIT/X style license                                   %
%    See License_ROMS.txt                           Hernan G. Arango        %
%===========================================================================%

% Set miscelaneous parameters.

global IPRINT
IPRINT=0;

format long g

% Change "my_root" to the appropriate path were these matlab scripts are
% installed in your computer.

my_home = getenv('HOME');
my_root = strcat(my_home, '/ocean/repository');

path(path, fullfile(my_root, 'matlab', '4dvar', ''))
path(path, fullfile(my_root, 'matlab', 'bathymetry', ''))
path(path, fullfile(my_root, 'matlab', 'boundary', ''))
path(path, fullfile(my_root, 'matlab', 'coastlines', ''))
path(path, fullfile(my_root, 'matlab', 'colormaps', ''))
path(path, fullfile(my_root, 'matlab', 'forcing', ''))
path(path, fullfile(my_root, 'matlab', 'grid', ''))
path(path, fullfile(my_root, 'matlab', 'grid_gui', ''))
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

% Load NetCDF Toolbox for OpenDAP support for versions 2008b or higher. 
% However, this is not needed if version 2012a or higher since Matlab
% native NetCDF interface supports OpenDAP.  Users need to change the
% paths for SNCTOOLS and JAVA.

v = version('-release');
vyear = str2num(v(1:4));
load_toolbox = vyear >= 2008;
if ((vyear == 2008 && v(5:5) == 'a') || vyear >= 2012),
  load_toolbox = false;
end

if (load_toolbox),
  addpath (strcat(my_home, '/ocean/matlab/snctools'), '-end');
  javaaddpath (strcat(my_home, '/ocean/matlab/classes/toolsUI-4.1.jar'), '-end');
  javaaddpath (strcat(my_home, '/ocean/matlab/classes/netcdfAll-4.2.jar'), '-end');
  javaaddpath (strcat(my_home, '/ocean/matlab/snctools/classes'), '-end');
  setpref('SNCTOOLS','USE_JAVA', true);
end