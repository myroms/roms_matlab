%
% ROMS Miscellaneous Matlab scripts
% =================================
%
% This utility contains several generic Matlab scripts to pre- and
% post-processing ROMS data.
%
% The NetCDF scripts use the MEXNC toolbox that can be downloaded
% from:
%
%   svn checkout url ~/matlab/mexnc
%
% where url = https://mexcdf.svn.sourceforge.net/svnroot/mexcdf/mexnc/trunk
%
% Similarly, the SNCTOOLS uses the NetCDF Java interface for OpenDAP
% that can be downloaded from:
%
%   svn checkout url ~/matlab/snctools
%
% where url = https://mexcdf.svn.sourceforge.net/svnroot/mexcdf/snctools/trunk
%
% Terrain-Following Coordinates:
%
%   depths        - Computes ROMS depths associated with a NetCDF 3D variable.
%   hslope        - Computes and plot ROMS grid bathymetry slope.
%   rfactor       - Computes bathymetry stiffness ratio, r-factor.
%   scoord        - Computes and plot ROMS vertical stretched coordinates.
%   set_depth     - Computes ROMS depths for 3D variable during pre-processing,
%                    like initial conditions, climatology, etc.
%   smooth_bath   - Smooths bathymetry as function of the r-factor.
%   stretching    - Computes ROMS vertical coordinate stretching function.
%
% NetCDF I/O Processing:
%
%   nc_attadd     - Adds/modifies a global or variable NetCDF attribute.
%   nc_attdel     - Deletes requested global or variable NetCDF attribute.
%   nc_dim        - Inquires about the dimensions in a NetCDF file.
%   nc_drename    - Renames a NetCDF dimension.
%   nc_getatt     - Gets a global or variable NetCDF attribute.
%   nc_inq        - Inquires about the contents of a NetCDF file.
%   nc_url        - Detects if the NetCDF file is a URL (OpenDAP server file).
%   nc_varinfo    - Returns structure with information about a NetCDF variable.
%   nc_vdef       - Creates a ROMS variable in a NetCDF file.
%   nc_vinfo      - Inquires information about requested NetCDF variable.
%   nc_vname      - Gets the names of all variables in a NetCDF file.
%   nc_vrename    - Renames a NetCDF variable.
%
%   nc_read       - Generic function to read  requested NetCDF variable.
%   nc_write      - Generic function to write requested NetCDF variable.
%
%   nc_slice      - Interpolates requested slice from a 3D NetCDF variable.
%
% ROMS Data Processing:
%
%   roms_vectors  - Processes vector data for either the full grid
%                     or boundary edges. The strategy is to get any
%                     horizontal vector field at RHO-points for the
%                     event that a rotation to ROMS curvilinear grid
%                     is needed.  Then, they are computed at the
%                     appropriate Arakawa C-grid location.
%
%   uv_barotropic - Computes vertically integrated velocity components
%                      for ROMS full grid or boundaries.
%
% Filters:
%
%   shapiro1      - 1D Shapiro filter.
%   shapiro2      - 2D Shapiro filter.
%
% Geophysical:
%
%   eos           - ROMS equation of state for seawater.
%   gcircle       - Great circle distance between two (lon,lat) points.
%   geodesic_dist - Geodesic distance between two (lon,lat) points.
%
% Time Management:
%
%   caldate      - Converts Julian day number to calendar date structure.
%   date_stamp   - Sets current date string.
%   day_code     - Computes day of the week for a given date.
%   gregorian    - Converts Julian day number to Gregorian calendar date.
%   greg2str     - Converts Gregorian date array to string.
%   hms2h        - Converts hours, minutes, and seconds to decimal hours.
%   julian       - Converts Gregorian calendar date to Julian day numbers.
%   s2hms        - Converts decimal seconds to integer hour, minute, seconds.
%
% Parallelism:
%
%   ptile        - Plot (overlay) ROMS horizontal tile partitions.
%   tile         - Compute ROMS horizontal tile partitions indices.
%

% svn $Id$
%===========================================================================%
%  Copyright (c) 2002-2011 The ROMS/TOMS Group                              %
%    Licensed under a MIT/X style license                                   %
%    See License_ROMS.txt                           Hernan G. Arango        %
%===========================================================================%
