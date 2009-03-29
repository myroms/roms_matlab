%
% ROMS Miscellaneous Matlab scripts
% =================================
%
% This utility contains several generic Matlab scripts to pre- and
% post-processing ROMS data.
%
%
% Terrain-Following Coordinates:
%
%   depths        - Compute ROMS depths associated with a NetCDF 3D variable.
%   scoord        - Compute and plot ROMS vertical stretched coordinates.
%   set_depth     - Compute ROMS depths for 3D variable during pre-processing,
%                    like initial conditions, climatology, etc.
%   stretching    - Compute ROMS vertical coordinate stretching function.
%
% NetCDF I/O Processing:
%
%   nc_dim        - Inquire about the dimensions in a NetCDF file.
%   nc_inq        - Inquire about the contents of a NetCDF file.
%   nc_vdef       - Create a ROMS variable in a NetCDF file.
%   nc_vinfo      - Inquire information about requested NetCDF variable.
%   nc_vname      - Get the names of all variables in a NetCDF file.
%
%   nc_read       - Generic function to read  requested NetCDF variable.
%   nc_write      - Generic function to write requested NetCDF variable.
%
%   nc_slice      - Interpolate requested slice from a 3D NetCDF variable.
%
% Preprocessing:
%
%   c_climatology - Create a ROMS climatology NetCDF file.
%   d_climatology - Driver to create a ROMS climatology file.
%
%   c_initial     - Create a ROMS initial conditions NetCDF file.
%   d_initial     - Driver to create a ROMS initial conditions file.
%
% Time management:
%
%   date_stamp   - Set current date string.
%   day_code     - Compute day of the week for a given date.
%
% Parallelism:
%
%   ptile        - Plot (overlay) ROMS horizontal tile partitions.
%   tile         - Compute ROMS horizontal tile partitions indices.
%

% svn $Id$
%===========================================================================%
%  Copyright (c) 2002-2009 The ROMS/TOMS Group                              %
%    Licensed under a MIT/X style license                                   %
%    See License_ROMS.txt                           Hernan G. Arango        %
%===========================================================================%
