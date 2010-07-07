%
% ROMS 4D-Var data assimilation Matlab scripts
% ===========================================
%
% This package contains several generic Matlab scripts to process
% data for ROMS 4D-Var data assimilation algorithms.
%
% 4D-Var Observations:
%
%   c_observations - Creates 4D-Var observation NetCDF file.
%   d_observations - Driver to process 4D-Var observation NetCDF file.
%
%   obs_read       - Reads observation NetCDF file and load all data
%                      into a structure array.
%   obs_write      - Writes all observation data in structure array into
%                      an existing NetCDF file.
%
%   inside         - Checks if a point is strictly inside of the region
%                      defined by a polygon. This is an old Matlab
%                      function which is not longer supported. It is
%                      very useful to find outliers observations.
%   obs_ijpos      - Computes observations locations in ROMS fractional
%                      coordinates. It uses the 'inside' deprecated
%                      Matlab function.
%
%   super_obs      - Checks the provided observation data (4D-Var NetCDF
%                      file or structure) and creates super observations if
%                      there is more than one datum of the same observation
%                      type per grid cell.
%   plot_super     - Sample script to compute and plot super observations
%                      from specified 4D-Var observations NetCDF file.
%
%   d_ssh_obs      - Driver template to extract SSH observations for
%                      AVISO, create and write observation file. Then,
%                      it computes super observations and creates and
%                      write super observations NetCDF file.
%   load_ssh_data  - Extracts AVISO sea level anomaly for the period of
%                      interest and specified region from ROMS Grid file.
%
% Error Covariance Matrix:
%
%   average:       - Computes the time average of requested NetCDF variable.
%   variance:      - Computes the variance of requested NetCDF variable from
%                     its specified time mean.
%
% Error Covariance Matrix Balance Operator:
%
%   balance_4dvar  - Computes 4D-Var balance operator.
%   biconj         - Biconjugate gradient solver for the SSH elliptic
%                      equation.
%   ini_balance    - Initializes balance operator structure array.
%                      It sets internal parameters, reads needed grid
%                      metrics and computes several quantities.
%   lateral_obc    - Sets lateral boundary conditions for a 2D or 3D
%                      field.
%   rho_balance    - Computes balanced density anomaly using a linear
%                      equation of state.
%   s_balance      - Given a temperature anomaly, deltaT=T-Tavg, it 
%                      computes balanced salinity anomaly using a T-S
%                      empirical formula.
%   ssh_reference  - Computes the balance operator reference sea surface
%                      height.
%   uv_balance     - Computes balanced, baroclinic U- and V-momentum
%                      anomalies (m/s) using the geostrophic balance.
%   zeta_balance   - Computes balanced, baroclinic free-surface anomaly
%                      by solving an elliptical equation OR integrating
%                      the hydrostatic equation from surface to bottom.
%
% Error Covariance Matrix NetCDF file:
%
%   c_std         - Creates 4D-Var standard deviation NetCDF file.
%

% svn $Id$
%===========================================================================%
%  Copyright (c) 2002-2010 The ROMS/TOMS Group                              %
%    Licensed under a MIT/X style license                                   %
%    See License_ROMS.txt                           Hernan G. Arango        %
%===========================================================================%
