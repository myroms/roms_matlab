function startup

% startup -- User script configuration for Matlab.  It can set default
%            paths, define Handle Graphics defaults, or predefine
%            variables in your workspace.
%
% NOTES:
%
% * The 'seagrid' GUI is deprecated. It only works for Matlab version
%   R2017a or lower.
%
% * The OpenDAP support for NetCDF files residing in remote web servers
%   is available for Matlab version R2012a or higher in the native NetCDF
%   interface.
%
% * The IODA conversion script from ROMS native 4D-Var to JEDI's IODA
%   format only works for Matlab version R2022b or higher because it
%   uses variable type NC_STRING (NetCDF4 files only).
%
% * We highly recommend that Users define the ROMS_ROOT_DIR variable in
%   their computer shell logging environment, specifying where the User
%   cloned/downloaded the ROMS source code, Test Cases, and Matlab
%   processing software:
%
%   setenv ROMS_ROOT_DIR  MyDownlodLocationDirectory

% svn $Id$
%=========================================================================%
%  Copyright (c) 2002-2023 The ROMS/TOMS Group                            %
%    Licensed under a MIT/X style license                                 %
%    See License_ROMS.txt                           Hernan G. Arango      %
%=========================================================================%

% Set miscelaneous parameters.

global IPRINT
IPRINT=0;

format long g

v = version('-release');
vyear = str2num(v(1:4));

% Change "my_root" to the appropriate path were these matlab scripts are
% installed in your computer.

my_root = getenv('ROMS_ROOT_DIR');

path(path, fullfile(my_root, 'roms_matlab', '4dvar', ''))
path(path, fullfile(my_root, 'roms_matlab', 'bathymetry', ''))
path(path, fullfile(my_root, 'roms_matlab', 'boundary', ''))
path(path, fullfile(my_root, 'roms_matlab', 'coastlines', ''))
path(path, fullfile(my_root, 'roms_matlab', 'colormaps', ''))
path(path, fullfile(my_root, 'roms_matlab', 'coupling', ''))
path(path, fullfile(my_root, 'roms_matlab', 'forcing', ''))
path(path, fullfile(my_root, 'roms_matlab', 'grid', ''))
path(path, fullfile(my_root, 'roms_matlab', 'initial', ''))
path(path, fullfile(my_root, 'roms_matlab', 'ioda', ''))
path(path, fullfile(my_root, 'roms_matlab', 'landmask', ''))
path(path, fullfile(my_root, 'roms_matlab', 'm_map', ''))
path(path, fullfile(my_root, 'roms_matlab', 'mex', ''))
path(path, fullfile(my_root, 'roms_matlab', 'netcdf', ''))
path(path, fullfile(my_root, 'roms_matlab', 'seagrid', ''))
path(path, fullfile(my_root, 'roms_matlab', 'seagrid', 'presto', ''))
path(path, fullfile(my_root, 'roms_matlab', 'seawater', ''))
path(path, fullfile(my_root, 'roms_matlab', 't_tide', ''))
path(path, fullfile(my_root, 'roms_matlab', 'tidal_ellipse', ''))
path(path, fullfile(my_root, 'roms_matlab', 'utility', ''))
