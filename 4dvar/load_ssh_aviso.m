function [D]=load_ssh_aviso(GRDfile, StartDay, EndDay)

%
% LOAD_SSH_AVISO:  Loads AVISO data for the region and time period
%
% [D]=load_ssh_aviso(GRDfile, StartDay, EndDay)
% 
%  Given a ROMS grid NetCDF, this function loads the AVISO sea level
%  anomaly covering the region for the time period specified (StartDay
%  to EndDay).
%   
% On Input:
%
%    GRDname       NetCDF grid file name (string)
%
%    StartDay      Starting period of interest (date number)
%
%                    Example:   StartDay=datenum(2004,1,1)=731947
%
%    EndDay        Ending   period of interest (date number)
%
%                    Example:   StartDay=datenum(2004,2,1)=731978
%
%                    See Matlab intrinsic datenum(Y,Mo,D,H,Mi,S)
%                    for details.
%
% On Output:
%
%    D            AVISO sea level data (structure array):
%
%                    D.ssh      sea level anomalies, cm (time,lat,lon)
%                    D.time     time of analysis (date number)
%                    D.lon      longitude of extracted data (lat,lon)
%                    D.lat      latitude  of extracted data (lat,lon)
%
% Warning: This function uses 'nc_varget' from SNCTOOLS with OpenDAP
%          to read NetCDF data.
%

% svn $Id$
%===========================================================================%
%  Copyright (c) 2002-2013 The ROMS/TOMS Group                              %
%    Licensed under a MIT/X style license                                   %
%    See License_ROMS.txt                           Brian Powell            %
%===========================================================================%

%  Check arguments.

if (nargin < 3),
  error([' LOAD_SSH_DATA: You must specify a grid file along with', ...
         ' starting and ending times']);
end,

if (StartDay > EndDay),
  error([' LOAD_SSH_DATA: Your starting time must be greater than', ...
         ' the ending time']);
end,

data=[];

%----------------------------------------------------------------------------
%  Extract AVISO data for the period of interest.
%----------------------------------------------------------------------------

%  Source of AVISO data.

if (StartDay < datenum(2010,3,31)),
  url=['http://opendap.aviso.oceanobs.com/thredds/dodsC/' ...
       'dataset-duacs-dt-ref-global-merged-msla-h?%s'];
else,
  url=['http://opendap.aviso.oceanobs.com/thredds/dodsC/' ...
       'duacs_global_nrt_msla_merged_h?%s'];
end,

file=sprintf(url,'time,NbLatitudes,NbLongitudes,Grid_0001');

%  Find the time period of interest.

epoch = datenum(1950,1,1);
ssh_time = nc_varget(file,'time') + epoch;
times = find(ssh_time >= StartDay & ssh_time <= EndDay);

if (isempty(times)),
  disp([' LOAD_SSH_DATA: no data found for period of interest.']);
  return;
end

ssh_time = ssh_time(times);
times = [min(times)-1 length(times)];

%  Find the region of interest.

ssh_lat = nc_varget(file,'NbLatitudes');
ssh_lon = nc_varget(file,'NbLongitudes');

%  Read in application grid longitude and latitude.

rlon = nc_varget(GRDfile,'lon_rho');
rlat = nc_varget(GRDfile,'lat_rho');

MinLon = min(rlon(:))-0.5;
MaxLon = max(rlon(:))+0.5;
MinLat = min(rlat(:))-0.5;
MaxLat = max(rlat(:))+0.5;

%  Check how western longitudes are handled.

ind = find(ssh_lon > MaxLon);
if (~isempty(ind)),
  ssh_lon(ind) = ssh_lon(ind) - 360;
end,

ind = find(ssh_lon < MinLon);
if (~isempty(ind)),
  ssh_lon(ind) = ssh_lon(ind) + 360;
end,

% Grab the indices for the application grid.

lon_list = find( ssh_lon >= MinLon & ssh_lon <= MaxLon );
lat_list = find( ssh_lat >= MinLat & ssh_lat <= MaxLat );

if (isempty(lon_list) | isempty(lat_list))
  return;
end,

[ssh_lon, ssh_lat] = meshgrid(ssh_lon(lon_list), ssh_lat(lat_list));

lon_list = [min(lon_list)-1 length(lon_list)];
lat_list = [min(lat_list)-1 length(lat_list)];


% Now that we have everything, load the data

ssh = nc_varget(file, 'Grid_0001', [times(1) lon_list(1) lat_list(1)], ...
                                   [times(2) lon_list(2) lat_list(2)]);

ssh = permute(ssh,[1 3 2]);        % ssh(time,lat,lon);

% We now have the times, position, and data. We'll put it into a structure
% for further processing

D.ssh  = ssh;
D.time = ssh_time;
D.lon  = ssh_lon;
D.lat  = ssh_lat;

return
