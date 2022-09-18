function [S]=ioda_read(ncfile)

%
% IODA_READ:  Reads IODA observation NetCDF4 file
%
% [S]=ioda_read(ncfile)
%
% This function reads IODA observation NetCDF4 file and stores all
% the variables in structure array, S.
%
% On Input:
%
%    ncfile  IODA Observations NetCDF4 file name (string)
%
% On Output:
%
%    S       IODA Observations data (structure array):
%
%              S.ncfile           NetCDF file name (string)
%              S.ndatetime        number of datetime
%              S.nlocs            number of locations
%              S.nobs             number of observations
%              S.nrecs            number of records
%              S.nstring          number of string
%              S.nvars            number of variables
%              S.datetime_ref     IODA reference time YYYYMMDDHH
%              S.variable_names   UFO/IODA standard name
%              S.dateTime         seconds since yyyy-mm-ddTHH:MM:SSZ
%              S.date_time        date and time ISO 8601 UTC string
%              S.datenum          date number
%              S.depth            depth of observations
%              S.longitude        longitude of observations
%              S.latitude         latitude  of observations
%              S.provenance       observation origin identifier
%              S.values           observation values
%              S.units            observation units
%              S.errors           observation error
%              S.PreQC            observation quality control
%

% svn $Id$
%=========================================================================%
%  Copyright (c) 2002-2022 The ROMS/TOMS Group                            %
%    Licensed under a MIT/X style license                                 %
%    See License_ROMS.txt                           Hernan G. Arango      %
%=========================================================================%
  
% Initialize.

Obs = struct('ncfile'        , [],                                      ...
             'ndatetime'     , [],                                      ...
             'nlocs'         , [],                                      ...
             'nobs'          , [],                                      ...
             'nrecs'         , [],                                      ...
             'nstring'       , [],                                      ...
             'nvars'         , [],                                      ...
             'datetime_ref'  , [],                                      ...
             'source'        , [],                                      ...
             'variable_names', [],                                      ...
             'date_time'     , [],                                      ...
             'depth'         , [],                                      ...
             'latitude'      , [],                                      ...
             'longitude'     , [],                                      ...
             'provenance'    , [],                                      ...
             'value'         , [],                                      ...
             'units'         , [],                                      ...
             'error'         , []);

% Inquire NetCDF4 file.

I = ncinfo(ncfile);

S.ncfile = ncfile;
S.ndatetime  = I.Dimensions(strcmp({I.Dimensions.Name}, 'ndatetime')).Length;
S.nlocs      = I.Dimensions(strcmp({I.Dimensions.Name}, 'nlocs'    )).Length;
S.nobs       = I.Dimensions(strcmp({I.Dimensions.Name}, 'nobs'     )).Length;
S.nrecs      = I.Dimensions(strcmp({I.Dimensions.Name}, 'nrecs'    )).Length;
S.nstring    = I.Dimensions(strcmp({I.Dimensions.Name}, 'nstring'  )).Length;
S.nvars      = I.Dimensions(strcmp({I.Dimensions.Name}, 'nvars'    )).Length;

% Get IODA reference time YYYYMMDDHH global attribute.

S.datetime_ref  = I.Attributes(strcmp({I.Attributes.Name}, 'date_time'  )).Value;
S.source        = I.Attributes(strcmp({I.Attributes.Name}, 'source_file')).Value;

% Read in 'VarMetaData' Group.

Vnames = transpose(ncread(ncfile, '/VarMetaData/variable_names'));

S.variable_names = {Vnames};

% Read in 'MetaData' Group.

G = ncinfo(ncfile, 'MetaData');

S.longitude = double(ncread(ncfile, '/MetaData/longitude'));
S.latitude  = double(ncread(ncfile, '/MetaData/latitude'));

if (any(strcmp({G.Variables.Name}, 'depth')))
  S.depth   = double(ncread(ncfile, '/MetaData/depth'));
end

S.date_time = transpose(ncread(ncfile, '/MetaData/date_time'));
S.datenum   = datenum(S.date_time, 'yyyy-mm-ddTHH:MM:SSZ');

S.timeDate = double(ncread(ncfile, '/MetaData/time'));

% Read in 'ObsValue' Group.

for i = 1:S.nvars
  Vname = ['/ObsValue/' S.variable_names{i}];
  values = double(ncread(ncfile, Vname));
  S.values(i) = {values};
  S.units(i) = {nc_getatt(ncfile, 'units', Vname)};
end

% Read in 'ObsError' Group.

for i = 1:S.nvars
  Vname = ['/ObsError/' S.variable_names{i}];
  errors = double(ncread(ncfile, Vname));
  S.errors(i) = {errors};
end

% Read in 'ObsError' Group.

for i = 1:S.nvars
  Vname = ['/PreQC/' S.variable_names{i}];
  preqc = double(ncread(ncfile, Vname));
  S.PreQC(i) = {preqc};
end

return
