%
% ROMS Initial Conditions and Climatology
% =======================================
%
% These functions are used for preparing ROMS initial conditions
% and climatology NetCDF file.
%
%
%   c_biology       - Defines ROMS biology variables to an existing
%                       initial conditions or climatology NetCDF file.
%   c_climatology   - Creates ROMS climatology NetCDF file.
%   c_initial       - Creates ROMS initial conditions NetCDF file. 
%
%   d_climatology   - Driver to process generic ROMS climatology file.
%   d_initial       - Driver ro process generic ROMS initial conditions
%                       file.
%   d_mercator2roms - Driver template to process ROMS initial conditions
%                       using the 'mercator2roms' script.
%   d_roms2roms     - Driver template to process ROMS initial conditions
%                       using the 'roms2roms' script.
%
%   mercator2roms - Interpolates requested variable data from
%                     Mercator to ROMS grids.
%   roms2roms     - Interpolates requested variable data from
%                     ROMS to ROMS grids.

% svn $Id$
%=========================================================================%
%  Copyright (c) 2002-2012 The ROMS/TOMS Group                            %
%    Licensed under a MIT/X style license                                 %
%    See License_ROMS.txt                           Hernan G. Arango      %
%=========================================================================%
