%
% ROMS Miscellaneous Matlab scripts
% ================================
%
% This utility contains several generic Matlab scripts to pre- and
% post-processing ROMS data.
%
%
% Terrain-Following Coordinates:
%
%   depths       - Compute ROMS depths associated with a 3D variable.
%   scoord       - Compute and plot ROMS vertical stretched coordinates.
%
% NetCDF I/O Processing:
%
%   nc_dim       - Inquire about the dimensions in a NetCDF file.
%   nc_inq       - Inquire about the contents of a NetCDF file.
%   nc_vdef      - Create a ROMS variable in a NetCDF file.
%   nc_vinfo     - Inquire information about requested NetCDF variable.
%   nc_vname     - Get the names of all variables in a NetCDF file.
%
%   nc_read      - Generic function to read  requested NetCDF variable.
%   nc_write     - Generic function to write requested NetCDF variable.
%
%   nc_slice     - Interpolate requested slice from a 3D NetCDF variable
%

% svn $Id$
%===========================================================================%
%  Copyright (c) 2002-2009 The ROMS/TOMS Group                              %
%    Licensed under a MIT/X style license                                   %
%    See License_ROMS.txt                           Hernan G. Arango        %
%===========================================================================%
