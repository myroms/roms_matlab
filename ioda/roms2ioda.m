function roms2ioda(ObsData, HisName, prefix, suffix, varargin)

%
% ROMS2IODA:  Converts ROMS observation NetCDF to IODA NetCDF4 files
%
% roms2ioda(ObsData, HisName, prefix, suffix, ...
%           SSHareaAvg, SSHtimeAvg, UVtimeAvg)
%
% This function converts a ROMS 4D-Var observation NetCDF file into several
% IODA NetCDF files. One file per each observation type is usually the way
% that the JEDI/UFO observation operator requires. Output files are of the
% form:
%                        prefix_obstype_suffix.nc4
%
% For example:           wc13_sst_20040103.nc4
%
% The area-averaged and time-averaged values will be used in ROMS H(x)
% computations with area-averaged and/or time-averaged filters.
%
% On Input:
%
%    ObsData     ROMS standard 4D-Var observation NetCDF filename (string)
%             or ROMS observation data structure (struct)
%
%    HisName     ROMS application history NetCDF filename (string)
%
%    prefix      Output file prefix associated with application (string)
%
%    suffix      Output file suffix associated with date and time, use
%                  YYYYMMDD or YYYYMMDDhh (numeric or string).
%
%                  NOTICE that the time of the output observations are
%                  converted to 'seconds since' time of suffix.
%
%                  For example, if suffix =  20040103  or
%                                  suffix = '20040103' then
%
%                  int64 dateTime(Location);
%                  dateTime:units = "seconds since 2004-01-03T00:00:00Z";
%
%                  Therefore, input and output time of the observations
%                  may have different time references!
%
%    SSHareaAvg  SSH observations area-averaged radius (km, OPTIONAL)
%                  Use unique SSH observations (remove repetitive)
%
%                  float spatialAverage
%
%    SSHtimeAvg  SSH observations time-averaged window since reference
%                  time (hours; OPTIONAL).
%
%                  For example 24 hours, 2004-01-03 to 2004-01-04
%                                        2004-01-04 to 2004-01-05
%                                        2004-01-05 to 2004-01-06
%                                        and so on ...
%
%                  int64 dateTimeAverageStart(timeWindow);
%
%                  int64 dateTimeAverageEnd(timeWindow);
%
%    UVtimeAvg   CODAR observations time-averaged window since reference
%                  time (hours; OPTIONAL)
%
%                  For example 24 hours, 2004-01-03 to 2004-01-04
%                                        2004-01-04 to 2004-01-05
%                                        2004-01-05 to 2004-01-06
%                                        and so on ...
%
%                  int64 dateTimeAverageStart(timeWindow);
%
%                  int64 dateTimeAverageEnd(timeWindow);
%
% Usage: Example to convert WC13 4D-Var observation file to create the
%        following IODA-2 files:
%                                  wc13_adt_20040103.nc4
%                                  wc13_sst_20040103.nc4
%                                  wc13_salt_20040103.nc4
%                                  wc13_temp_20040103.nc4
%
%   ObsData = 'WC13/Data/wc13_obs_20040103.nc';
%   HisName = 'WC13/Forward/r06/wc13_roms_his_20040103.nc';
%
%   roms2ioda(ObsData, HisName, 'wc13', '20040103')
%
% Dependencies:
%
%                create_ioda_obs.m
%                get_roms_grid.m
%                obs_k2z.m
%                obs_read.m

% svn $Id$
%=========================================================================%
%  Copyright (c) 2002-2025 The ROMS Group                                 %
%    Licensed under a MIT/X style license                                 %
%    See License_ROMS.md                            Hernan G. Arango      %
%=========================================================================%

% Initialize.

switch numel(varargin)
  case 0
    SSHareaAvg = NaN;
    SSHtimeAvg = NaN;
    UVtimeAvg  = NaN;
  case 1
    SSHareaAvg = varargin{1};
    SSHtimeAvg = NaN;
    UVtimeAvg  = NaN;
  case 2
    SSHareaAvg = varargin{1};
    SSHtimeAvg = varargin{2};
    UVtimeAvg  = NaN;
  case 3
    SSHareaAvg = varargin{1};
    SSHtimeAvg = varargin{2};
    UVtimeAvg  = varargin{3};
end

if (ischar(ObsData))
  S = obs_read(ObsData);
else
  S = ObsData;
end
I = ncinfo(S.ncfile);              % NetCDF information structure

% Get ROMS grid structure.

G = get_roms_grid(HisName, HisName);

% Build IODA date and time numeric (YYYYMMDDhh) global attribute
% 'date_time'.  It is used as reference starting date and time to
% compute the MetaData Group 'time' coordinate in hours.

if (isnumeric(suffix))
  suffix = num2str(suffix);
end

if (length(suffix) == 8)
  DateTimeIODA = str2num([suffix '00']);
else
  DateTimeIODA = str2num(suffix);
end

% Check if the 'reference_time' field in observation structure is
% available and not empty.

if (~isfield(S, 'reference_time'))
  error(['ROMS2IODA: cannot find ''reference_time'' in observation'     ...
         ' structure. Get update ''obs_read.m'' function']);
end

if (isempty(isfield(S, 'reference_time')))
  error(['ROMS2IODA: cannot find ''reference_date'' in time unit'       ...
         ' attribute. Non compliant file: ', S.ncfile])
end

% Make sure that observation depths are in negative meters.

S.Zgrid = S.depth;
S = obs_k2z(S, HisName);

% Estimate data assimilation cycle time-window.

days_window = floor((max(S.time)-min(S.time))+0.5);

%--------------------------------------------------------------------------
% Inquire observation structure about type meassurent types.
%--------------------------------------------------------------------------

% Get switches for observations for state variables.

types = unique(S.type);

got_ssh  = any(types == 1);
got_uvel = any(types == 4);
got_vvel = any(types == 5);
got_temp = any(types == 6);
got_salt = any(types == 7);

% Get indices for state variables observations.

issh  = find(S.type == 1);
iuvel = find(S.type == 4);
ivvel = find(S.type == 5);
itemp = find(S.type == 6);
isalt = find(S.type == 7);

% Currently, one of the strategies to improve the impact of altimetry
% track data is to repeat their observations several times over the
% assimilation cycle, but with different errors.
%
% Identify such repetitive observations by setting their provenance
% to negative values.

if (~isempty(issh))
  [~,IA,~] = unique(complex(S.lat(issh), S.lon(issh)));
  if (~isempty(IA))
    sshProv = -abs(S.provenance(issh));      % Set negative SSH provenance
    for n = 1:length(IA)
      sshProv(IA(n))= -sshProv(IA(n));       % Turn positive unique values
    end
    S.provenance(issh) = sshProv;            % overwrite SSH provenance
  end

  if (~isnan(SSHareaAvg) || ~isnan(SShtimeAbg))
    issh = find(S.provenance(issh) > 0);
  end
end

% Get provenance indices for each state variable.

P = struct('ssh', [], 'uvel', [], 'vvel', [], 'temp', [], 'salt', []);

if (~isempty(issh))
  P.ssh = unique(S.provenance(issh));
end
if (~isempty(iuvel))
  P.uvel = unique(S.provenance(iuvel));
end
if (~isempty(ivvel))
  P.vvel = unique(S.provenance(ivvel));
end
if (~isempty(itemp))
  P.temp = unique(S.provenance(itemp));
end
if (~isempty(isalt))
  P.salt = unique(S.provenance(isalt));
end

% If available, get the 'flag_values' and 'flag_meanings' attributes
% for 'obs_provenance' variable.

got_flagAtt = false;
if (any(strcmp({I.Variables(strcmp({I.Variables.Name},                  ...
               'obs_provenance')).Attributes.Name}, 'flag_values')))
  flag_values = nc_getatt(S.ncfile, 'flag_values', 'obs_provenance');
  flag_meanings = nc_getatt(S.ncfile, 'flag_meanings', 'obs_provenance');

  P.flag_values   = flag_values;
  P.flag_meanings = strsplit(flag_meanings);
  got_flagAtt = true;
end

% Process sea surface height observations. IODA wants absolute
% dynamic topography (ADT).
%
% ADT = MDT + SLA   (Mean Dynamic Topography plus Sea Level Anomaly)
% SSH = MSS + SLA   (Mean Sea Surface plus Sea Level Anomaly)

if (got_ssh)
  if (~isempty(issh))
    has_depth  = false;
    Obs = extract_observations(S, issh, DateTimeIODA, has_depth);
    Obs.ncfile         = [prefix '_adt_' suffix '.nc4'];
    Obs.N              = G.N;
    Obs.nvars          = 1;
    Obs.units          = {'meter'};
    Obs.ncvname        = {'absoluteDynamicTopography'};
    Obs.stateID        = 1;
    Obs.areaAvgRadius  = SSHareaAvg;
    Obs.timeAvgWindow  = SSHtimeAvg;
    Obs.variables_name = {'absolute_dynamic_topography'};
    Obs.datetime_ref   = DateTimeIODA;
    if (got_flagAtt)
      Obs.flag_values   = int32(P.flag_values(P.ssh));
      Obs.flag_meanings = string(join(P.flag_meanings(P.ssh)));
    end
    if (~isnan(SSHareaAvg))
      Obs.areaAvgRadius = SSHareaAvg * 1000;        % to meters
    end
    if (~isnan(SSHtimeAvg))
      delta  = SSHtimeAvg * 3600;                    % to secods
      window = days_window * 86400;
      Tstr = 0:delta:window-delta;
      Tend = delta:delta:window;
      Obs.timeAvgBegin = Tstr;
      Obs.timeAvgEnd   = Tend;
      Obs.nwindow      = length(Tstr);
    end
    create_ioda_obs(Obs);
    write_observations(Obs);
  end
end

% Process sea surface temperature.

if (got_temp)
  isst = find(S.type == 6 & S.Zgrid == G.N);
  if (~isempty(isst))
    has_depth = false;
    Obs = extract_observations(S, isst, DateTimeIODA, has_depth);
    Obs.ncfile         = [prefix '_sst_' suffix '.nc4'];
    Obs.N              = G.N;
    Obs.nvars          = 1;
    Obs.units          = {'C'};
    Obs.ncvname        = {'seaSurfaceTemperature'};
    Obs.stateID        = 6;
    Obs.variables_name = {'sea_surface_temperature'};
    Obs.datetime_ref   = DateTimeIODA;
    if (got_flagAtt)
      sst_ind = contains(upper(P.flag_meanings), 'SST');
      if (any(sst_ind))
        Obs.flag_values   = int32(P.flag_values(sst_ind));
        Obs.flag_meanings = string(join(P.flag_meanings(sst_ind)));
      else
        Obs.flag_values   = int32(P.flag_values(P.temp));
        Obs.flag_meanings = string(join(P.flag_meanings(P.temp)));
      end
    end
    create_ioda_obs(Obs);
    write_observations(Obs);
  end
end

% Process sub-surface insitu temperature.

if (got_temp)
  ktemp = find(S.type == 6 & S.Zgrid ~= G.N);
  if (~isempty(ktemp))
    has_depth = true;
    Obs = extract_observations(S, ktemp, DateTimeIODA, has_depth);
    Obs.ncfile         = [prefix '_temp_' suffix '.nc4'];
    Obs.N              = G.N;
    Obs.nvars          = 1;
    Obs.units          = {'C'};
    Obs.stateID        = 6;
    Obs.ncvname        = {'waterTemperature'};
    Obs.variables_name = {'sea_water_temperature'};
    Obs.datetime_ref   = DateTimeIODA;
    if (got_flagAtt)
      sst_ind = contains(upper(P.flag_meanings(P.temp)), 'SST');
      if (any(sst_ind))                             % remove SST provenance
        P.temp(sst_ind)   = [];
        Obs.flag_values   = int32(P.flag_values(P.temp));
        Obs.flag_meanings = string(join(P.flag_meanings(P.temp)));
      else
        Obs.flag_values   = int32(P.flag_values(P.temp));
        Obs.flag_meanings = string(join(P.flag_meanings(P.temp)));
      end
    end
    create_ioda_obs(Obs);
    write_observations(Obs);
  end
end

% Process sub-surface potential temperature.  It is assumed that all the
% observations in the input NetCDF are insitu temperature. Below, they
% are converted to potential temperature.

if (got_temp)
  ktemp = find(S.type == 6 & S.Zgrid ~= G.N);
  if (~isempty(ktemp))
    has_depth = true;
    Obs = extract_observations(S, ktemp, DateTimeIODA, has_depth);

    p = sw_pres(abs(Obs.depth), Obs.latitude);      % decibars
    s = ones(size(Obs.ObsValue)).*35.0;
    pr = zeros(size(Obs.ObsValue));
    ptemp = sw_ptmp(s, Obs.ObsValue, p, pr);
    Obs.ObsValue = ptemp;

    Obs.ncfile         = [prefix '_ptemp_' suffix '.nc4'];
    Obs.N              = G.N;
    Obs.nvars          = 1;
    Obs.units          = {'C'};
    Obs.ncvname        = {'waterPotentialTemperature'};
    Obs.stateID        = 6;
    Obs.variables_name = {'sea_water_potential_temperature'};
    Obs.datetime_ref   = DateTimeIODA;
    if (got_flagAtt)
      sst_ind = contains(upper(P.flag_meanings(P.temp)), 'SST');
      if (any(sst_ind))                             % remove SST provenance
        P.temp(sst_ind)   = [];
        Obs.flag_values   = int32(P.flag_values(P.temp));
        Obs.flag_meanings = string(join(P.flag_meanings(P.temp)));
      else
        Obs.flag_values   = int32(P.flag_values(P.temp));
        Obs.flag_meanings = string(join(P.flag_meanings(P.temp)));
      end
    end
    create_ioda_obs(Obs);
    write_observations(Obs);
  end
end

% Process salinity.

if (got_salt)
  if (~isempty(isalt))
    has_depth = true;
    Obs = extract_observations(S, isalt, DateTimeIODA, has_depth);
    Obs.ncfile         = [prefix '_salt_' suffix '.nc4'];
    Obs.N              = G.N;
    Obs.nvars          = 1;
    Obs.units          = {'dimensionless'};
    Obs.ncvname        = {'salinity'};
    Obs.stateID        = 7;
    Obs.variables_name = {'sea_water_salinity'};
    Obs.datetime_ref   = DateTimeIODA;
    if (got_flagAtt)
      Obs.flag_values   = int32(P.flag_values(P.salt));
      Obs.flag_meanings = string(join(P.flag_meanings(P.salt)));
    end
    create_ioda_obs(Obs);
    write_observations(Obs);
  end
end

% Process gridded HF radar near surface (z=-2 m) velocities. The
% CODAR velocity components in the native ROMS observation file
% are rotated to grid curvilinear coordinates, we need to rotate
% back to geographycal EAST and NORTH directions before they are
% written to the IODA-type file.

if (got_uvel && got_vvel)
  if ~(isempty(iuvel) || isempty(iuvel))
    has_depth = true;
    Uobs = extract_observations(S, iuvel, DateTimeIODA, has_depth);
    Vobs = extract_observations(S, ivvel, DateTimeIODA, has_depth);
    Obs  = Uobs;
    Obs.ObsError       = [];
    Obs.ObsError{1}    = Uobs.ObsError;
    Obs.ObsError{2}    = Vobs.ObsError;
    Obs.ObsValue       = [];
    Obs.ObsValue{1}    = Uobs.ObsValue;
    Obs.ObsValue{2}    = Vobs.ObsValue;
    Obs.PreQC          = [];
    Obs.PreQC{1}       = Uobs.PreQC;
    Obs.PreQC{2}       = Uobs.PreQC;
    Obs.ncfile         = [prefix '_uv_codar_' suffix '.nc4'];
    Obs.N              = G.N;
    Obs.nvars          = 2;
    Obs.units          = {'m s-1', 'm s-1'};
    Obs.ncvname        = {'waterZonalVelocity',                      ...
                          'waterMeridionalVelocity'};
    Obs.stateID        = [4, 5];
    Obs.timeAvgWindow  = UVtimeAvg;
    Obs.variables_name = {'eastward_sea_water_velocity',             ...
                          'meridional_sea_water_velocity'};
    Obs.datetime_ref   = DateTimeIODA;
    if (got_flagAtt)
      Obs.flag_values   = int32(P.flag_values(P.uvel));
      Obs.flag_meanings = string(join(P.flag_meanings(P.uvel)));
    end
    if (~isnan(UVtimeAvg))
      delta  = UVtimeAvg * 3600;                    % to secods
      window = days_window * 86400;
      Tstr = 0:delta:window-delta;
      Tend = delta:delta:window;
      Obs.timeAvgBegin = Tstr;
      Obs.timeAvgEnd   = Tend;
      Obs.nwindow      = length(Tstr);
    end
    Obs = rotate_codar(Obs, G);
    create_ioda_obs(Obs);
    write_observations(Obs);
  end
end

return

%--------------------------------------------------------------------------
% Extract requested observation type from structure.

function [Obs] = extract_observations(S, ind, DateTimeIODA, has_depth)

% Initialize structure.

Obs = struct('ncfile'        , [],                                      ...
             'source'        , [],                                      ...
             'N'             , [],                                      ...
             'nlocs'         , [],                                      ...
             'nobs'          , [],                                      ...
             'nsurvey'       , [],                                      ...
             'nvars'         , [],                                      ...
             'nwindow'       , NaN,                                     ...
             'units'         , [],                                      ...
             'ncvname'       , [],                                      ...
             'variable_names', [],                                      ...
             'TimeIODA'      , [],                                      ...
             'DateIODA'      , [],                                      ...
             'dateTime'      , [],                                      ...
             'surveyTime'    , [],                                      ...
             'surveyIndex'   , [],                                      ...
             'date_time'     , [],                                      ...
             'depth'         , [],                                      ...
             'latitude'      , [],                                      ...
             'longitude'     , [],                                      ...
             'provenance'    , [],                                      ...
             'stateID'       , [],                                      ...
             'sequenceNumber', [],                                      ...
             'areaAvgRadius' , NaN,                                     ...
             'timeAvgBegin'  , NaN,                                     ...
             'timeAvgEnd'    , NaN,                                     ...
             'x_grid'        , [],                                      ...
             'y_grid'        , [],                                      ...
             'z_grid'        , [],                                      ...
             'ObsError'      , [],                                      ...
             'ObsValue'      , [],                                      ...
             'PreQC'         , []);

% Fill some values in the structure.

Obs.source         = S.ncfile;

Obs.nlocs          = length(ind);
Obs.nvars          = 1;
Obs.nobs           = Obs.nlocs * Obs.nvars;
if (has_depth)
  Obs.depth        = S.depth(ind);
end
Obs.latitude       = S.lat(ind);
Obs.longitude      = S.lon(ind);
Obs.provenance     = S.provenance(ind);
Obs.sequenceNumber = int32(1:Obs.nlocs);
Obs.x_grid         = S.Xgrid(ind);
Obs.y_grid         = S.Ygrid(ind);
if (has_depth)
  Obs.z_grid       = S.Zgrid(ind);
end
Obs.ObsError       = sqrt(S.error(ind));              % standard deviation
Obs.ObsValue       = S.value(ind);
Obs.PreQC          = int32(zeros(size(Obs.longitude)));

% Compute dateTime seconds from IODA file reference time.

ioda_epoch      = datenum(num2str(DateTimeIODA), 'yyyymmddHH');

time_days       = S.reference_time + S.time(ind);
ioda_time_days  = time_days - ioda_epoch;

Obs.dateTime    = int64(ioda_time_days * 86400);      % seconds

[survey,IA,IC]  = unique(Obs.dateTime, 'stable');
Obs.surveyTime  = int64(survey);
Obs.surveyIndex = int32(IA);
Obs.nsurvey     = length(Obs.surveyTime);

Obs.date_time   = cellstr(datestr(ioda_epoch + ioda_time_days,          ...
                                 'yyyy-mm-ddTHH:MM:SSZ'));

% Set IODA file time reference global attributes.

Obs.TimeIODA = DateTimeIODA;
Obs.DateIODA = datestr(ioda_epoch, 'yyyy-mm-ddTHH:MM:SSZ');

return

%--------------------------------------------------------------------------
% Extract requested observation type from structure.

function [Obs] = rotate_codar (S, G)

% Interpolate ROMS curvilinear angle to CODAR observation locations.
% Use scatteredInterpolant since horizontal grid coordinates are not
% plaid.

F = scatteredInterpolant(G.lon_rho(:), G.lat_rho(:), G.angle(:));

angle = F(S.longitude, S.latitude);

% Rotate velocity component to geographical EAST and NORTH.

Urot = S.ObsValue{1};
Vrot = S.ObsValue{2};

cos_angle = cos(angle);
sin_angle = sin(angle);

Ugeo = Urot .* cos_angle - Vrot .* sin_angle;
Vgeo = Vrot .* cos_angle + Urot .* sin_angle;

% Replace velocity components.

Obs = S;

Obs.ObsValue{1} = Ugeo;
Obs.ObsValue{2} = Vgeo;

return

%--------------------------------------------------------------------------
% Writes requested observation type data to ouput NetCDF file.

function write_observations(Obs)

disp(['*** Writing  observations file:  ', Obs.ncfile]);

% Write out 'MetaData' Group variables.

ncwrite(Obs.ncfile, 'MetaData/dateTime', int64(Obs.dateTime));
if (~isnan(Obs.timeAvgBegin))
  ncwrite(Obs.ncfile, 'MetaData/dateTimeAverageBegin',                  ...
          int64(Obs.timeAvgBegin));
end
if (~isnan(Obs.timeAvgEnd))
  ncwrite(Obs.ncfile, 'MetaData/dateTimeAverageEnd',                    ...
          int64(Obs.timeAvgEnd));
end
if (~isempty(Obs.depth))
  ncwrite(Obs.ncfile, 'MetaData/depth', Obs.depth);
end
ncwrite(Obs.ncfile, 'MetaData/latitude', Obs.latitude);
ncwrite(Obs.ncfile, 'MetaData/longitude', Obs.longitude);
if (~isempty(Obs.provenance))
  ncwrite(Obs.ncfile, 'MetaData/provenance', Obs.provenance);
end
ncwrite(Obs.ncfile, 'MetaData/stateID', int32(Obs.stateID));
ncwrite(Obs.ncfile, 'MetaData/sequenceNumber', Obs.sequenceNumber);
if (~isnan(Obs.areaAvgRadius))
  ncwrite(Obs.ncfile, 'MetaData/spatialAverage', Obs.areaAvgRadius);
end
ncwrite(Obs.ncfile, 'MetaData/surveyIndex', int32(Obs.surveyIndex));
ncwrite(Obs.ncfile, 'MetaData/surveyTime', int64(Obs.surveyTime));
ncwrite(Obs.ncfile, 'MetaData/variables_name', Obs.variables_name);
ncwrite(Obs.ncfile, 'MetaData/x_grid', Obs.x_grid);
ncwrite(Obs.ncfile, 'MetaData/y_grid', Obs.y_grid);
if (~isempty(Obs.z_grid))
  ncwrite(Obs.ncfile, 'MetaData/z_grid', Obs.z_grid);
end

% Write out 'ObsError' Group variables.

if (iscell(Obs.ObsError))
  for i = 1:Obs.nvars
    Vname = ['/ObsError/' Obs.ncvname{i}];
    ncwrite(Obs.ncfile, Vname, Obs.ObsError{i});
  end
else
  Vname = ['/ObsError/' char(Obs.ncvname)];
  ncwrite(Obs.ncfile, Vname, Obs.ObsError);
end

% Write out 'ObsValue' Group variables.

if (iscell(Obs.ObsValue))
  for i = 1:Obs.nvars
    Vname = ['/ObsValue/' Obs.ncvname{i}];
    ncwrite(Obs.ncfile, Vname, Obs.ObsValue{i});
  end
else
  Vname = ['/ObsValue/' char(Obs.ncvname)];
  ncwrite(Obs.ncfile, Vname, Obs.ObsValue);
end

% Write out 'ObsQC' Group variables.

if (iscell(Obs.PreQC))
  for i = 1:Obs.nvars
    Vname = ['/PreQC/' Obs.ncvname{i}];
    ncwrite(Obs.ncfile, Vname, Obs.PreQC{i});
  end
else
  Vname = ['/PreQC/' char(Obs.ncvname)];
  ncwrite(Obs.ncfile, Vname, Obs.PreQC);
end

return
