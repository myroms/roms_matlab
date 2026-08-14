function ioda_add_filter(InpName, OutName, M);

%
% IODA_ADD_FILTER:  Adds area-averaging or time-averaging filter
%
% This function duplicates input IODA NetCDF-4 observation file and
% variables in the 'MetaData' Group for area- and time-averaging 
% processing by the H(x) operators in a ROMS data assimilation
% application.
%
% int64 dateTimeAverageBegin(Location) ;
%  	dateTimeAverageBegin:long_name = "start of time averaging filter" ;
%  	dateTimeAverageBegin:units = "seconds since YYYY-MM-DDThh:mm:ssZ" ;
%	dateTimeAverageBegin:filter = "??-hour half-length averaging" ;
%
% int64 dateTimeAverageEnd(Location) ;
%       dateTimeAverageEnd:long_name = "end of time averaging filter" ;
%       dateTimeAverageEnd:units = "seconds since YYYY-MM-DDThh:mm:ssZ" ;
%       dateTimeAverageEnd:filter = "??-hour half-length averaging" ;
%
% float spatialAverage(nvars) ;
%       spatialAverage:long_name = "spatial averaging radius" ;
%       spatialAverage:units = "meter" ;  
%
% On Input:
%
%    InpName     Input IODA NetCDF-4 filename (string)
%
%    OutName     Output thinned IODA NetCDF-4 filename (string)
%
%    M           IODA enhanced NetCDF-4 file Metadata structure
%                  (struct array) computed with "ioda_metadat.m"
%
%                  M(:).name           variable short name
%                  M(:).cycle_length   Data Assimilation cyle (hours)
%                  M(:).radius         area-averaged radius (km) scale
%                  M(:).time_window    time-averaged window (hours)
%                  M(:).ioda_vname     IODA NetCDF-4 variable name
%                  M(:).standard_name  variable standard name
%
% USAGE:
% *****
%
%   Example to append area and time averaged filter to SSH obesrvations:
%
%   (1) Set IODA NetCDF-4 file metadata structure, M:
%
%       M = ioda_metadata(true);
%
%   (2) Set Data Assimilation cycle length (hours) using "deal" to
%       assing values to all structure elements.
%
%       [M.cycle_length] = deal(96);                             % hours
%
%   (3) If applicable, set area-averaged and time-averaged parameters
%       for specialized H(x) operators.
%
%       M(strcmp({M.name}, 'SSH')).radius = 30;                  % km
%       M(strcmp({M.name}, 'SSH')).time_window = 12;             % hours
%
%       To display updated values, use:
%
%         disp( struct2table(M) );
%
% git $Id$
%=========================================================================%
%  Copyright (c) 2002-2026 The ROMS Group                                 %
%    Licensed under a MIT/X style license                                 %
%    See License_ROMS.md                            Hernan G. Arango      %
%=========================================================================%

% Initialize area-averaged and time-averaged parameters from Metdata
% structure, M.

SSHareaAvg = M(strcmp({M.name}, 'SSH')).radius;
SSHtimeAvg = M(strcmp({M.name}, 'SSH')).time_window;

SSTareaAvg = M(strcmp({M.name}, 'SST')).radius;
SSTtimeAvg = M(strcmp({M.name}, 'SST')).time_window;

SSSareaAvg = M(strcmp({M.name}, 'SSS')).radius;
SSStimeAvg = M(strcmp({M.name}, 'SSS')).time_window;

UVareaAvg  = M(strcmp({M.name}, 'uv_CODAR')).radius;
UVtimeAvg  = M(strcmp({M.name}, 'uv_CODAR')).time_window;

SALTareaAvg = M(strcmp({M.name}, 'salt')).radius;
SALTtimeAvg = M(strcmp({M.name}, 'salt')).time_window;

TEMPareaAvg = M(strcmp({M.name}, 'temp')).radius;
TEMPtimeAvg = M(strcmp({M.name}, 'temp')).time_window;

% Read native ROMS observation file and load to S structure.

if (ischar(InpName))
  S = ioda_read(InpName);
else
  S = InpName;
end

% Get data assimilation cycle time-window.

if ~isempty([M.cycle_length])
  days_window = unique([M.cycle_length]) / 24;
else
  days_window = floor((max(S.dateTime)/86400)+0.5);
end

%------------------------------------------------------------------------
%  Append filter metadata to IODA observation structure.
%------------------------------------------------------------------------

switch S.ncvname{1}
  case 'absoluteDynamicTopography'
    if (~isempty(SSHtimeAvg))
      S.areaAvgRadius = SSHareaAvg * 1000;          % to meters
    end
    if (~isempty(SSHtimeAvg))
      delta  = SSHtimeAvg * 3600;                   % to secods
      window = days_window * 86400;
      Tstr = min(max(0, S.dateTime-delta), window);
      Tend = min(max(0, S.dateTime+delta), window);
      S.timeAvgBegin = Tstr;
      S.timeAvgEnd   = Tend;
      S.average_window = SSHtimeAvg;
    end
  case 'seaSurfaceTemperature'
    if (~isempty(SSTareaAvg))
      S.areaAvgRadius = SSTareaAvg * 1000;          % to meters
    end
    if (~isempty(SSTtimeAvg))
      delta  = SSTtimeAvg * 3600;                   % to secods
      window = days_window * 86400;
      Tstr = min(max(0, S.dateTime-delta), window);
      Tend = min(max(0, S.dateTime+delta), window);
      S.timeAvgBegin = Tstr;
      S.timeAvgEnd   = Tend;
      S.average_window = SSTtimeAvg;
    end 
  case 'seaSurfaceSalinity'
    if (~isempty(SSSareaAvg))
      S.areaAvgRadius = SSSareaAvg * 1000;           % to meters
    end
    if (~isempty(SSStimeAvg))
      delta  = SSStimeAvg * 3600;                    % to secods
      window = days_window * 86400;
      Tstr = min(max(0, S.dateTime-delta), window);
      Tend = min(max(0, S.dateTime+delta), window);
      S.timeAvgBegin = Tstr;
      S.timeAvgEnd   = Tend;
      S.average_window = SSStimeAvg;
    end 
  case {'waterZonalVelocity', 'waterSurfaceMeridionalVelocity'}
    if (~isempty(UVareaAvg))
      S.areaAvgRadius = UVareaAvg * 1000;           % to meters
    end 
    if (~isempty(UVtimeAvg))
      delta  = UVtimeAvg * 3600;                    % to secods
      window = days_window * 86400;
      Tstr = min(max(0, S.dateTime-delta), window);
      Tend = min(max(0, S.dateTime+delta), window);
      Obs.timeAvgBegin = Tstr;
      Obs.timeAvgEnd   = Tend;
      Obs.average_window = UVtimeAvg;
    end
  case {'waterTemperature', 'waterPotentialTemperature'}
    if (~isempty(TEMPtimeAvg))
      delta  = TEMPtimeAvg * 3600;                   % to secods
      window = days_window * 86400;
      Tstr = min(max(0, S.dateTime-delta), window);
      Tend = min(max(0, S.dateTime+delta), window);
      S.timeAvgBegin = Tstr;
      S.timeAvgEnd   = Tend;
      S.average_window = TEMPtimeAvg;
    end 
  case 'salinity'
    if (~isempty(SALTtimeAvg))
      delta  = SALTtimeAvg * 3600;                   % to secods
      window = days_window * 86400;
      Tstr = min(max(0, S.dateTime-delta), window);
      Tend = min(max(0, S.dateTime+delta), window);
      S.timeAvgBegin = Tstr;
      S.timeAvgEnd   = Tend;
      S.average_window = SALTtimeAvg;
    end 
  otherwise
    error(['Cannot process variable: ', S.ncvname{1}]);
end

%  Create IODA NetCDF-4 file.

S.ncfile = OutName;
create_ioda_obs(S);

%  Write data into IODA NetCDF-4 file.

ioda_write(S, OutName);

return

