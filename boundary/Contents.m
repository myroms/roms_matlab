%
% ROMS Open Boundary Conditions
% =============================
%
% These functions are used for preparing ROMS open boundary conditions
% NetCDF file.
%
%
%   c_boundary    - Creates ROMS open boundary conditions NetCDF
%                     file.
%
%   d_mercator    - Driver to extract open boundary conditions from
%                     Mercator dataset. It creates open boundary
%                     conditions NetCDF file, interpolates data to
%                     application grid, and writes out data.
%
%   obc_mercator  - Interpolates requested variable open boundary
%                     conditions from Mercator to ROMS grid boundary
%                     edges.
%
%   extract_bry   - Reads requested variable from a ROMS NetCDF file
%                     at the specified time record and extracts the
%                     lateral boundary edges. No interpolation is
%                     carried out.
%

% svn $Id$
%===========================================================================%
%  Copyright (c) 2002-2011 The ROMS/TOMS Group                              %
%    Licensed under a MIT/X style license                                   %
%    See License_ROMS.txt                           Hernan G. Arango        %
%===========================================================================%
