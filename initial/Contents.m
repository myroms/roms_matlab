%
% ROMS Initial Conditions and Climatology
% =======================================
%
% These functions are used for preparing ROMS initial conditions
% and climatology NetCDF file.
%
%
%   c_biology     - Defines ROMS biology variables to an existing
%                     initial conditions or climatology NetCDF file.
%   c_climatology - Creates ROMS climatology NetCDF file.
%   c_initial     - Creates ROMS initial conditions NetCDF file. 
%
%   d_climatology - Driver to process ROMS climatology file.
%   d_initial     - Driver ro process ROMS initial conditions file.
%
%   mercator2roms - Interpolates requested variable data from
%                     Mercator to ROMS grid.
%

% svn $Id$
%===========================================================================%
%  Copyright (c) 2002-2010 The ROMS/TOMS Group                              %
%    Licensed under a MIT/X style license                                   %
%    See License_ROMS.txt                           Hernan G. Arango        %
%===========================================================================%