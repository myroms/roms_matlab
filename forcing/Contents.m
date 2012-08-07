%
% ROMS Forcing Fields
% ===================
%
% These functions are used for preparing ROMS forcing NetCDF files.
%
%
%   d_core2_frc   - Driver template script showing how to create ROMS
%                     forcing NetCDF file(s) using ROMS metadata
%                     structure. The data source is the CORE 2 Global
%                     Air-Sea Flux Dataset. Notice that the original
%                     data set is sampled for the Gulf of Mexico (GOM).
%
%  write_tides    - Creates ROMS tidal forcing NetCDF file and writes
%                     data extracted from either OTPS or ADCIRC and
%                     processed with the "t_tide" utility.
%

% svn $Id$
%=========================================================================%
%  Copyright (c) 2002-2012 The ROMS/TOMS Group                            %
%    Licensed under a MIT/X style license                                 %
%    See License_ROMS.txt                           Hernan G. Arango      %
%=========================================================================%