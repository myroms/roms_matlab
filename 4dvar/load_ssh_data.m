function [data]=load_ssh_data(GRDfile,startday,endday)

%
% LOAD_SSH_DATA:  Loads AVISO data for the region and time period
%
% [DATA]=load_ssh_data(GRDfile,startday,endday)
% 
%  Given a ROMS grid NetCDF, this function loads the AVISO sea level
%  anomaly covering the region for the time period specified (startday
%  to endday).
%   
% On Input:
%
%    GRDname       NetCDF grid file name (string)
%
%    startday      starting period of interest (date number)
%
%                    Example:   startday=datenum(2004,1,1)=731947
%
%    startday      ending   period of interest (date number)
%
%                    Example:   startday=datenum(2004,2,1)=731978
%
%                    See Matlab intrinsic datenum(Y,Mo,D,H,Mi,S)
%                    for details.
%
% On Output:
%
%    data          AVISO sea level anomalies (structure array):
%
%                    Data.ssh      sea level anomalies
%                    Data.time     time of analysis (date number)
%                    Data.lon      longitude of extracted data
%                    Data.lat      latitude  of extracted data
%
% Warning: This function uses 'nc_varget' from SNCTOOLS with OpenDAP
%          to read NetCDF data.
%

% svn $Id$
%===========================================================================%
%  Copyright (c) 2002-2010 The ROMS/TOMS Group                              %
%    Licensed under a MIT/X style license                                   %
%    See License_ROMS.txt                           Brian Powell            %
%===========================================================================%

%  Check arguments.

if (nargin < 3),
  error([' LOAD_SSH_DATA: You must specify a grid file along with', ...
	 ' starting and ending times']);
end,

if (startday > endday),
  error([' LOAD_SSH_DATA: Your starting time must be greater than', ...
	 ' the ending time']);
end,

data=[];

%----------------------------------------------------------------------------
%  Extract AVISO data for the period of interest.
%----------------------------------------------------------------------------

%  Source of AVISO data.

if (startday < datenum(2001,8,22)),
  url='http://opendap.aviso.oceanobs.com/thredds/dodsC/dt_ref_global_merged_msla_h?%s';
else,
  url='http://opendap.aviso.oceanobs.com/thredds/dodsC/duacs_global_nrt_msla_merged_h?%s';
end,

file=sprintf(url,'time,NbLatitudes,NbLongitudes,Grid_0001');

%  Find the time period of interest.

epoch=datenum(1950,1,1);
ssh_time=nc_varget(file,'time') + epoch;
times=find(ssh_time >= startday & ssh_time <= endday);
if (isempty(times)),
  disp([' LOAD_SSH_DATA: no data found for period of interest.']);
  return;
end
ssh_time=ssh_time(times);
times=[min(times)-1 length(times)];

%  Find the region of interest.

ssh_lat=nc_varget(file,'NbLatitudes');
ssh_lon=nc_varget(file,'NbLongitudes');

%  Read in application grid longitude and latitude.
	
rlon=nc_varget(GRDfile,'lon_rho');
rlat=nc_varget(GRDfile,'lat_rho');
	
% Check how western longitudes are handled.

mnlon=min(rlon(:))-0.5;
mxlon=max(rlon(:))+0.5;
mnlat=min(rlat(:))-0.5;
mxlat=max(rlat(:))+0.5;

if (mxlon < 0),
  l=find(ssh_lon > 180);
  ssh_lon(l)=ssh_lon(l)-360;
end,

% Grab the indices for the application grid.

lon_list = find( ssh_lon >= mnlon & ssh_lon <= mxlon );
lat_list = find( ssh_lat >= mnlat & ssh_lat <= mxlat );

if (isempty(lon_list) | isempty(lat_list))
  return;
end,

[ssh_lon, ssh_lat]=meshgrid(ssh_lon(lon_list), ssh_lat(lat_list));

lon_list=[min(lon_list)-1 length(lon_list)];
lat_list=[min(lat_list)-1 length(lat_list)];


% Now that we have everything, load the data

ssh = nc_varget(file, 'Grid_0001', [times(1) lon_list(1) lat_list(1)], ...
                                   [times(2) lon_list(2) lat_list(2)]);

ssh = permute(ssh,[1 3 2]);        % ssh(time,lat,lon);

% We now have the times, position, and data. We'll put it into a structure
% for further processing

data.ssh  = ssh;
data.time = ssh_time;
data.lon  = ssh_lon;
data.lat  = ssh_lat;

return
