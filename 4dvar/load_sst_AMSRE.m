function [data]=load_sst_AMSRE(GRDfile, StartDay, EndDay)

%
% LOAD_SST_AMSRE:  Loads SST data for the region and time period
%
% [DATA]=load_sst_pfeg(GRDfile, StartDay, EndDay)
% 
%  Given a ROMS grid NetCDF, this function loads the SST data from the
%  extensive OpenDAP catalog maintained by NOAA PFEG OceanWatch. The
%  SST has 0.025 degree global 1-day product. Meassurements are gathered
%  by Japan's Advanced Microwave Scanning Radiometer (AMSR-E) instrument,
%  a passive radiance sensor carried aboard NASA's Aqua spacecraft.
%  The entire OceanWatch catalog can be found at:  
%  
%  http://oceanwatch.pfeg.noaa.gov/thredds/catalog.html
%
%  This script process the SST, Aqua AMSR-E, Global/1-day dataset which
%  starts on 2002-06-01 12:00:00Z
%
%  http://oceanwatch.pfeg.noaa.gov/thredds/Satellite/aggregsatAA/ssta/catalog.html
%   
% On Input:
%
%    GRDname       NetCDF grid file name (string)
%
%    StartDay      Starting period of interest (date number)
%
%                    Example:   startday=datenum(2004,1, 1)=731947
%
%    EndDay        Ending   period of interest (date number)
%
%                    Example:   startday=datenum(2004,1,15)=731961
%
%                    See Matlab intrinsic datenum(Y,Mo,D,H,Mi,S)
%                    for details.
%
% On Output:
%
%    data          Composite SST data (structure array):
%
%                    Data.time     time of extracted data (date number)
%                    Data.lon      longitude of extracted data
%                    Data.lat      latitude  of extracted data
%                    Data.sst      sea surface temperatures
%
% Warning: This function uses 'nc_varget' from SNCTOOLS with OpenDAP
%          to read NetCDF data.
%

% svn $Id$
%===========================================================================%
%  Copyright (c) 2002-2013 The ROMS/TOMS Group                              %
%    Licensed under a MIT/X style license                                   %
%    See License_ROMS.txt                           John Wilkin             %
%===========================================================================%

%  Check arguments.

if (nargin < 3),
  error([' LOAD_SST_AMSRE: You must specify a grid file along with', ...
         ' starting and ending times']);
end,

if (StartDay > EndDay),
  error([' LOAD_SST_AMSRE: Your starting time must be greater than', ...
         ' the ending time']);
end,

data=[];

%----------------------------------------------------------------------------
%  Extract SST data for the period of interest.
%----------------------------------------------------------------------------

%  Source of SST data.

url = ['http://thredds1.pfeg.noaa.gov:8080/thredds/dodsC/' ...
       'satellite/AA/ssta/1day'];

%  Find the time period of interest.

epoch = datenum([1970 1 1 0 0 0]);
sst_time = epoch + nc_varget(url,'time')/86400;
T = find(sst_time >= StartDay & sst_time <= EndDay);

if (isempty(T)),
  disp([' LOAD_SST_AMSRE: no data found for period of interest.']);
  return;
end

T = T - 1;            %  substract 1 because SNCTOOLS is 0-based.

%  Read SST longitudes and latitudes.

sst_lon = nc_varget(url,'lon');
sst_lat = nc_varget(url,'lat');

%  Read in application grid longitude and latitude.

rlon = nc_varget(GRDfile,'lon_rho');
rlat = nc_varget(GRDfile,'lat_rho');

MinLon = min(rlon(:))-0.5;
MaxLon = max(rlon(:))+0.5;
MinLat = min(rlat(:))-0.5;
MaxLat = max(rlat(:))+0.5;

%  Check how western longitudes are handled.

ind = find(sst_lon > MaxLon);
if (~isempty(ind)),
  sst_lon(ind) = sst_lon(ind) - 360;
end,

ind = find(sst_lon < MinLon);
if (~isempty(ind)),
  sst_lon(ind) = sst_lon(ind) + 360;
end,

%  Grab the indices for the application grid.

I = find(sst_lon >= MinLon & sst_lon <= MaxLon);
J = find(sst_lat >= MinLat & sst_lat <= MaxLat);

if (isempty(I) | isempty(J))
  disp([' LOAD_SST_AMSRE: no data found for application grid.']);
  return;
end,

I = I - 1;            %  substract 1 because SNCTOOLS is 0-based.
J = J - 1;            %  substract 1 because SNCTOOLS is 0-based.

%  Read again coordinates for the selected region and time period to
%  be safe.

data.time = nc_varget(url,'time',T(1), length(T));
data.lon  = nc_varget(url,'lon' ,I(1), length(I));
data.lat  = nc_varget(url,'lat' ,J(1), length(J));

data.time = epoch + data.time./86400;

ind = find(data.lon > MaxLon);
if (~isempty(ind)),
  data.lon(ind) = data.lon(ind) - 360;
end,

ind = find(data.lon < MinLon);
if (~isempty(ind)),
  data.lon(ind) = data.lon(ind) + 360;
end,

%  Get the SST data. The data are actually 4D with second coordinate
%  being altitude.

data.sst = nc_varget(url, 'AAssta', [T(1) 0 J(1) I(1)], ...
                                    [length(T) 1 length(J) length(I)]);   

return
