%
%  PLOT_SUPER:  Plot super observations
%
%  This script can be used to check the effect of binning 4D-Var
%  observations after processing with 'super_obs'.
%

% svn $Id$
%===========================================================================%
%  Copyright (c) 2002-2010 The ROMS/TOMS Group                              %
%    Licensed under a MIT/X style license                                   %
%    See License_ROMS.txt                           Hernan G. Arango        %
%===========================================================================%

%  Set input/output NetCDF files.

 my_root = '/home/arango/ocean/toms/repository/test';

 GRDname = fullfile(my_root, 'WC13/Data', 'wc13_grd.nc');
 OBSname = 'wc13_ssh_obs.nc';

%  Read in application grid.

rlon=nc_read(GRDname,'lon_rho');
rlat=nc_read(GRDname,'lat_rho');

got_coast=true;
try,
  clon=nc_read(GRDname,'lon_coast');
  clat=nc_read(GRDname,'lat_coast');
catch
  got_coast=false;
end,

%  Read in observations into structure.

[O]=obs_read(OBSname);

%  Create super observations, if necessary.

[N]=super_obs(O);

%  Plot original (O) and new (N) observations.

state_var=1;                          % Check altimetry (AVISO)
Nsurvey=length(O.survey_time);        % Number of surveys
survey_time=O.survey_time(1);         % select first survey

ind_o=find(O.type == state_var & O.time == survey_time);

if (~isempty(ind_o)),
  Olon=O.lon(ind_o);
  Olat=O.lat(ind_o);
  Oval=O.value(ind_o);
end,

ind_n=find(N.type == state_var & N.time == survey_time);

if (~isempty(ind_n)),
  Nlon=N.lon(ind_n);
  Nlat=N.lat(ind_n);
  Nval=N.value(ind_n);
end,

pcolor(rlon,rlat,ones(size(rlon)).*NaN);
hold on;
if (got_coast),
  plot(clon,clat,'k-');
end,
scatter(Olon,Olat,30,Oval,'ko');
scatter(Nlon,Nlat,30,Nval,'d','filled');
colorbar('fontsize',10,'fontweight','bold');
hold off;
title('Original (black circles), Super Obs (filled diamonds)', ...
      'fontsize',10,'fontweight','bold');



