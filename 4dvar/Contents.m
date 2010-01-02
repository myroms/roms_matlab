%
% ROMS 4DVar data assimilation Matlab scripts
% ===========================================
%
% This package contains several generic Matlab scripts to process
% data for ROMS 4DVar data assimilation algorithms.
%
%
% Error Covariance Matrix:
%
%   average:      - Computes the time average of requested NetCDF variable.
%   variance:     - Computes the variance of requested NetCDF variable from
%                     its specified time mean.
%
% Error Covaraince Matrix Balance Operator:
%
%   balance_4dvar - Computes 4DVar balance operator.
%   biconj        - Biconjugate gradient solver for the SSH elliptic
%                     equation.
%   ini_balance   - Initializes balance operator structure array.
%                     It sets internal parameters, reads needed grid
%                     metrics and computes several quantities.
%   lateral_obc   - Sets lateral boundary conditions for a 2D or 3D
%                     field.
%   rho_balance   - Computes balanced density anomaly using a linear
%                     equation of state.
%   s_balance     - Given a temperature anomaly, deltaT=T-Tavg, it 
%                     computes balanced salinity anomaly using a T-S
%                     empirical formula.
%   ssh_reference - Computes the balance operartor reference sea surface
%                     height.
%   uv_balance    - Computes balanced, baroclinic U- and V-momentum
%                     anomalies (m/s) using the geostrophic balance.
%   zeta_balance  - Computes balanced, baroclinic free-surface anomaly
%                     by solving an elliptical equation OR integrating
%                     the hydrostatic equation from surface to bottom.
%
% Error Covariance Matrix NetCDF file:
%
%   c_std         - Creates 4DVar standard deviation NetCDF file.
%

% svn $Id$
%===========================================================================%
%  Copyright (c) 2002-2010 The ROMS/TOMS Group                              %
%    Licensed under a MIT/X style license                                   %
%    See License_ROMS.txt                           Hernan G. Arango        %
%===========================================================================%
